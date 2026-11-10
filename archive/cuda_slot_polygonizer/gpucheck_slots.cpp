// gpucheck_slots.cpp -- ARCHIVED 2026-11-09 (plan D2). Not built, not included. The dense slot-grid mask polygonizer
// of cuda/gpucheck.cpp as of commit b6c2df2 (D1), the last commit that used it: its masks section (MaskView, has_edge,
// step_x/y, Slots, InitSlots, LabelRound, Pred, RankInit, RankRound, RingStats, ScatterCorners, jump_rounds,
// kChunkSlots, GatherRings, SetBase, check_mask_chunk) and its check_masks entry. Replaced by the K1-draft port
// (docs/implementation2.0/kernels/K1) in cuda/gpucheck.cpp. See README.md.

// ===================================================================== masks
// Slot of the unit edge leaving vertex (x, y) in direction d (0 east, 1 south, 2 west, 3 north), local to a mask:
// (y (W + 1) + x) 4 + d, which is also polygonize's edge key, so the minimum slot of a ring is its canonical start.
// The bit-quad at vertex (x, y) has TL = pixel (x-1, y-1) = 1, TR = (x, y-1) = 2, BL = (x-1, y) = 4, BR = (x, y) = 8,
// as in checker.bitquad. Foreground lies right of every edge (clockwise around parts, counter-clockwise in holes).
struct MaskView {
  const u8* bits;
  const i64 *start, *slot_start, *vrow_start;
  const i32 *W, *H;
  i64 n;

  HD bool px(i64 m, i32 x, i32 y) const {
    return x >= 0 && y >= 0 && x < W[m] && y < H[m] && bits[start[m] + i64(y) * W[m] + x] != 0;
  }
  HD int quad(i64 m, i32 x, i32 y) const {
    return px(m, x - 1, y - 1) | px(m, x, y - 1) << 1 | px(m, x - 1, y) << 2 | px(m, x, y) << 3;
  }
};

HD inline bool has_edge(int code, int d) {
  switch (d) {
    case 0: return (code & 8) && !(code & 2);  // top side of BR
    case 1: return (code & 4) && !(code & 8);  // right side of BL
    case 2: return (code & 1) && !(code & 4);  // bottom side of TL
    default: return (code & 2) && !(code & 1);  // left side of TR
  }
}

HD inline i32 step_x(int d) { return d == 0 ? 1 : (d == 2 ? -1 : 0); }
HD inline i32 step_y(int d) { return d == 1 ? 1 : (d == 3 ? -1 : 0); }
HD inline bool is_corner(const i32* pred, i32 e) { return (e & 3) != (pred[e] & 3); }

struct Slots {     // dense per-slot arrays of one chunk of masks, indexed slot_start[m] + local slot
  i32* succ;       // successor slot, -1 = no edge
  i32* pred;
  i32 *lab[2], *jmp[2], *w[2];  // ping-pong buffers of the pointer-jumping rounds
  i32* rcount;     // per ring (at its start slot): corners
  i32* rkey;       // per ring: first east edge in raster order (polygonize's ring order)
  i64* rarea;      // per ring: twice the signed area
  i64* rbase;      // per ring: output offset of its first corner
};

// The outer loops below walk one vertex row (mask m, y) per item; f(m, base, x, y, local slot e).
template <class F>
HD inline void for_row_slots(const MaskView& v, i64 i, const Slots& sl, F f) {
  const i64 m = segment(v.vrow_start, v.n, i);
  const i32 y = i32(i - v.vrow_start[m]), W = v.W[m];
  const i64 base = v.slot_start[m];
  for (i32 x = 0; x <= W; ++x)
    for (int d = 0; d < 4; ++d) {
      const i32 e = (y * (W + 1) + x) * 4 + d;
      if (sl.succ[base + e] >= 0) f(m, base, x, y, e);
    }
}

struct InitSlots {
  MaskView v;
  Slots sl;
  i64* quads;
  HD void operator()(i64 i) const {
    const i64 m = segment(v.vrow_start, v.n, i);
    const i32 y = i32(i - v.vrow_start[m]), W = v.W[m];
    const i64 base = v.slot_start[m];
    i64 n1 = 0, n2 = 0, nd = 0, n3 = 0, n4 = 0;
    for (i32 x = 0; x <= W; ++x) {
      const int code = v.quad(m, x, y);
      const int bits = (code & 1) + (code >> 1 & 1) + (code >> 2 & 1) + (code >> 3);
      if (code == 6 || code == 9) ++nd;
      else if (bits == 1) ++n1;
      else if (bits == 2) ++n2;
      else if (bits == 3) ++n3;
      else if (bits == 4) ++n4;
      for (int d = 0; d < 4; ++d) {
        const i32 e = (y * (W + 1) + x) * 4 + d;
        i32 next = -1;
        if (has_edge(code, d)) {  // successor: the left turn if it exists (saddles), else the only edge
          const i32 x1 = x + step_x(d), y1 = y + step_y(d);
          const int c1 = v.quad(m, x1, y1);
          const int turns[3] = {(d + 3) & 3, d, (d + 1) & 3};
          for (int t : turns)
            if (next < 0 && has_edge(c1, t)) next = (y1 * (W + 1) + x1) * 4 + t;
        }
        sl.succ[base + e] = next;
        sl.lab[0][base + e] = e;
        sl.jmp[0][base + e] = next;
        sl.rkey[base + e] = 0x7FFFFFFF;
      }
    }
    i64* q = quads + 5 * m;
    atomic_add(q + 0, n1), atomic_add(q + 1, n2), atomic_add(q + 2, nd), atomic_add(q + 3, n3), atomic_add(q + 4, n4);
  }
};

// Label round: lab = min(lab, lab[jmp]); jmp = jmp[jmp]. After k rounds lab covers 2^k slots of the cycle.
struct LabelRound {
  MaskView v;
  Slots sl;
  int in;
  HD void operator()(i64 i) const {
    const i32 *lab = sl.lab[in], *jmp = sl.jmp[in];
    i32 *lab2 = sl.lab[1 - in], *jmp2 = sl.jmp[1 - in];
    for_row_slots(v, i, sl, [&](i64, i64 base, i32, i32, i32 e) {
      const i32 j = jmp[base + e];
      lab2[base + e] = mn(lab[base + e], lab[base + j]);
      jmp2[base + e] = jmp[base + j];
    });
  }
};

struct Pred {
  MaskView v;
  Slots sl;
  HD void operator()(i64 i) const {
    for_row_slots(v, i, sl, [&](i64, i64 base, i32, i32, i32 e) { sl.pred[base + sl.succ[base + e]] = e; });
  }
};

// Ranking on pred, stopped at the ring's start: w(e) = corners from the start up to e's predecessor = e's corner index.
struct RankInit {
  MaskView v;
  Slots sl;
  const i32* lab;
  HD void operator()(i64 i) const {
    for_row_slots(v, i, sl, [&](i64, i64 base, i32, i32, i32 e) {
      const bool start = lab[base + e] == e;
      const i32 p = sl.pred[base + e];
      sl.jmp[0][base + e] = start ? e : p;
      sl.w[0][base + e] = start ? 0 : is_corner(sl.pred + base, p);
    });
  }
};

struct RankRound {
  MaskView v;
  Slots sl;
  int in;
  HD void operator()(i64 i) const {
    const i32 *w = sl.w[in], *jmp = sl.jmp[in];
    i32 *w2 = sl.w[1 - in], *jmp2 = sl.jmp[1 - in];
    for_row_slots(v, i, sl, [&](i64, i64 base, i32, i32, i32 e) {
      const i32 j = jmp[base + e];
      w2[base + e] = w[base + e] + w[base + j];
      jmp2[base + e] = jmp[base + j];
    });
  }
};

struct RingStats {
  MaskView v;
  Slots sl;
  const i32* lab;
  i32 *nrings, *rec;  // rec: (mask chunk index, start slot) of every ring, in any order
  HD void operator()(i64 i) const {
    for_row_slots(v, i, sl, [&](i64 m, i64 base, i32 x, i32 y, i32 e) {
      const i32 root = lab[base + e], d = e & 3, x1 = x + step_x(d), y1 = y + step_y(d);
      atomic_add(sl.rarea + base + root, i64(x) * y1 - i64(x1) * y);
      if (is_corner(sl.pred + base, e)) atomic_add(sl.rcount + base + root, 1);
      if (d == 0) atomic_min(sl.rkey + base + root, y * (v.W[m] + 1) + x);
      if (root == e) {
        const i32 k = atomic_inc(nrings);
        rec[2 * k] = i32(m), rec[2 * k + 1] = e;
      }
    });
  }
};

struct ScatterCorners {
  MaskView v;
  Slots sl;
  const i32 *lab, *w;
  i32* pts;
  HD void operator()(i64 i) const {
    for_row_slots(v, i, sl, [&](i64, i64 base, i32 x, i32 y, i32 e) {
      if (!is_corner(sl.pred + base, e)) return;
      const i64 k = sl.rbase[base + lab[base + e]] + w[base + e];
      pts[2 * k] = x, pts[2 * k + 1] = y;
    });
  }
};

int jump_rounds(i64 n) {  // polygonize._jumps: ceil(log2(n)) + 1 rounds cover any cycle of <= n edges
  int k = 0;
  while ((i64(1) << k) < n) ++k;
  return k + 1;
}

constexpr i64 kChunkSlots = i64(1) << 23;  // slots per batch of masks (~0.5 GB of scratch)

struct RingInfo {  // one ring of a mask chunk, gathered for the host's ordering
  i32 m, e, key, count;  // mask (chunk index), start slot, first east edge, corners
  i64 area2;             // twice the signed area
};

struct GatherRings {
  MaskView v;
  Slots sl;
  const i32* rec;
  RingInfo* info;
  HD void operator()(i64 k) const {
    const i32 m = rec[2 * k], e = rec[2 * k + 1];
    const i64 slot = v.slot_start[m] + e;
    info[k] = RingInfo{m, e, sl.rkey[slot], sl.rcount[slot], sl.rarea[slot]};
  }
};

struct SetBase {  // each ring's output offset, at its start slot
  MaskView v;
  Slots sl;
  const i32* rec;
  const i64* base;
  HD void operator()(i64 k) const { sl.rbase[v.slot_start[rec[2 * k]] + rec[2 * k + 1]] = base[k]; }
};

// Masks [m0, m1) on one lane: copy in, rings and corner ranks, then three short round trips to the host (the ring
// count; the rings and bit-quads; the corners after the host ordered the rings).
void check_mask_chunk(const Masks& in, i64 m0, i64 m1, MaskResult& out, Lane& L, Timer& t) {
  const i64 n = m1 - m0;
  std::vector<i64> slot_start(n + 1, 0), vrow_start(n + 1, 0);
  i64 max_slots = 1;
  for (i64 k = 0; k < n; ++k) {
    const i64 W = in.W[m0 + k], H = in.H[m0 + k];
    slot_start[k + 1] = slot_start[k] + 4 * (W + 1) * (H + 1);
    vrow_start[k + 1] = vrow_start[k] + H + 1;
    max_slots = std::max(max_slots, 4 * (W + 1) * (H + 1));
  }
  const i64 npx = in.start[m1] - in.start[m0], S = slot_start[n], rows = vrow_start[n];
  const i64 rmax = S / 8 + 1;  // <= 2 edges per vertex (4 slots), >= 4 edges per ring
  Layout lay;
  const size_t o_bits = lay.put<u8>(npx), o_start = lay.put<i64>(n + 1), o_slot = lay.put<i64>(n + 1),
               o_vrow = lay.put<i64>(n + 1), o_W = lay.put<i32>(n), o_H = lay.put<i32>(n);
  const size_t in_bytes = lay.size;
  size_t o_i32[9];  // succ, pred, lab 0/1, jmp 0/1, w 0/1, rkey
  for (size_t& o : o_i32) o = lay.put<i32>(S);
  const size_t o_rbase = lay.put<i64>(S), o_rec = lay.put<i32>(2 * rmax), o_info = lay.put<RingInfo>(rmax),
               o_base = lay.put<i64>(rmax), o_pts = lay.put<i32>(S);  // <= S / 2 corners
  const size_t o_rcount = lay.put<i32>(S), o_rarea = lay.put<i64>(S), o_quads = lay.put<i64>(5 * n),
               o_nrings = lay.put<i32>(1);  // [o_rcount, end): zeroed
  const size_t info_bytes = (sizeof(RingInfo) * size_t(rmax) + 255) / 256 * 256;
  u8* d = L.dev.get(lay.size);
  u8* h = L.in.get(std::max(in_bytes, sizeof(i64) * size_t(rmax)));
  u8* o = L.out.get(std::max(info_bytes + sizeof(i64) * size_t(5 * n), sizeof(i32) * size_t(S)));
  std::copy_n(in.bits + in.start[m0], npx, ptr<u8>(h, o_bits));
  for (i64 k = 0; k <= n; ++k) ptr<i64>(h, o_start)[k] = in.start[m0 + k] - in.start[m0];
  std::copy_n(slot_start.data(), n + 1, ptr<i64>(h, o_slot));
  std::copy_n(vrow_start.data(), n + 1, ptr<i64>(h, o_vrow));
  std::copy_n(in.W + m0, n, ptr<i32>(h, o_W));
  std::copy_n(in.H + m0, n, ptr<i32>(h, o_H));

  const Stream st = L.stream;
  const MaskView v{ptr<u8>(d, o_bits), ptr<i64>(d, o_start), ptr<i64>(d, o_slot), ptr<i64>(d, o_vrow), ptr<i32>(d, o_W),
                   ptr<i32>(d, o_H), n};
  auto slots = [&](int k) { return ptr<i32>(d, o_i32[k]); };
  const Slots sl{slots(0), slots(1), {slots(2), slots(3)}, {slots(4), slots(5)}, {slots(6), slots(7)},
                 ptr<i32>(d, o_rcount), slots(8), ptr<i64>(d, o_rarea), ptr<i64>(d, o_rbase)};
  i32* rec = ptr<i32>(d, o_rec);
  to_device(d, h, in_bytes, st);
  zero(d + o_rcount, lay.size - o_rcount, st);
  const int rounds = jump_rounds(max_slots);
  const i32 *lab = nullptr, *w = nullptr;
  t.stage("ccl", st, [&] {  // ring labels: pointer jumping to each cycle's minimum slot
    for_each(st, rows, InitSlots{v, sl, ptr<i64>(d, o_quads)});
    int in_buf = 0;
    for (int k = 0; k < rounds; ++k, in_buf ^= 1) for_each(st, rows, LabelRound{v, sl, in_buf});
    lab = sl.lab[in_buf];
  });
  t.stage("rings", st, [&] {  // corner ranks and the ring records
    for_each(st, rows, Pred{v, sl});
    for_each(st, rows, RankInit{v, sl, lab});
    int in_buf = 0;
    for (int k = 0; k < rounds; ++k, in_buf ^= 1) for_each(st, rows, RankRound{v, sl, in_buf});
    w = sl.w[in_buf];
    for_each(st, rows, RingStats{v, sl, lab, ptr<i32>(d, o_nrings), rec});
  });
  to_host(o, d + o_nrings, sizeof(i32), st);
  wait(st);
  const i32 R = *ptr<i32>(o, 0);
  if (R < 0 || R > rmax) throw std::logic_error("gpucheck: ring count out of its bound");

  t.stage("rings", st, [&] { for_each(st, R, GatherRings{v, sl, rec, ptr<RingInfo>(d, o_info)}); });
  to_host(o, d + o_info, sizeof(RingInfo) * size_t(R), st);
  to_host(o + info_bytes, d + o_quads, sizeof(i64) * size_t(5 * n), st);
  wait(st);
  // Order the rings as polygonize does and give each its output offset (host: rings are few).
  const RingInfo* info = ptr<RingInfo>(o, 0);
  std::vector<i32> order(R);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](i32 a, i32 b) {
    return info[a].m != info[b].m ? info[a].m < info[b].m : info[a].key < info[b].key;
  });
  const i64 p0 = out.pts_start.back();
  i64 total = 0;
  std::vector<i64> per_mask(n, 0);
  i64* base = ptr<i64>(h, 0);  // the copy in is done (wait above): its staging buffer is free
  for (i32 k : order) {
    base[k] = total;
    total += info[k].count;
    out.pts_start.push_back(p0 + total);
    out.ring_area.push_back(info[k].area2 / 2);
    per_mask[info[k].m] += 1;
  }
  for (i64 k = 0; k < n; ++k) out.ring_start.push_back(out.ring_start.back() + per_mask[k]);
  const i64* quads = ptr<i64>(o + info_bytes, 0);
  out.quads.insert(out.quads.end(), quads, quads + 5 * n);
  if (2 * total > S) throw std::logic_error("gpucheck: corner count out of its bound");

  to_device(d + o_base, base, sizeof(i64) * size_t(R), st);
  t.stage("trace", st, [&] {  // every corner to its output position
    for_each(st, R, SetBase{v, sl, rec, ptr<i64>(d, o_base)});
    for_each(st, rows, ScatterCorners{v, sl, lab, w, ptr<i32>(d, o_pts)});
  });
  to_host(o, d + o_pts, sizeof(i32) * size_t(2 * total), st);
  wait(st);
  out.pts.insert(out.pts.end(), ptr<i32>(o, 0), ptr<i32>(o, 0) + 2 * total);
}

// ---------------------------------------------------------------- the entry (gpucheck.cpp API section)
MaskResult check_masks(const Masks& in, Times* times) {
  GC_RANGE("check_masks");
  MaskResult out;
  out.ring_start = {0};
  out.pts_start = {0};
  if (in.n == 0) return out;
  Engine& e = engine();
  std::lock_guard<std::mutex> lock(e.mu);
  clear_error();
  Drain drain{e};
  Timer t{times, {}};
  e.start();
  for (i64 m0 = 0; m0 < in.n;) {  // chunks of at most kChunkSlots slots (a larger mask goes alone)
    i64 m1 = m0, slots = 0;
    while (m1 < in.n) {
      const i64 s = 4 * (i64(in.W[m1]) + 1) * (i64(in.H[m1]) + 1);
      if (m1 > m0 && slots + s > kChunkSlots) break;
      slots += s, ++m1;
    }
    check_mask_chunk(in, m0, m1, out, e.lanes[0], t);
    m0 = m1;
  }
  if (times) t.read();
  return out;
}
