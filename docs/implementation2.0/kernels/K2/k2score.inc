// k2score.inc -- K2: exact scoring and evaluation kernels, written as a fragment of gpucheck.cpp's unnamed namespace.
// It uses gpucheck's platform section (for_each, Stream, the atomic helpers, popc / ctz), its helpers (segment,
// floordiv, iou) and its polygon structures (Shapes, Rings, Edge, classify, EdgePairs, ShapeTotals). k2.cpp includes
// gpucheck.cpp and then this file inside the same unnamed namespace, so every functor below is a drop-in for it.
//
// Same rule as gpucheck: each algorithm is written once, as __host__ __device__ functors. for_each runs one item per
// thread; for_each_warp (new here) runs one item per warp of 32 lanes on the GPU and one item per "warp of 1 lane" on
// the CPU build, so the lane loops below (for i = lane; i < n; i += size) become plain loops and the warp collectives
// (wballot, wany, wodd, wsum, wmin) become the identity. Every output is an integer, a boolean, or a float64 made by
// ONE thread from integers in a fixed order (no FMA contraction: CMakeLists.txt), so any schedule, any atomic order and
// either build give the same bits.
//
//  1. ShapeTotalsWarp: gpucheck's ShapeTotals with one warp per shape. Ring sums and flags take one lane per ring; a
//     witness search ("does any vertex or edge midpoint of ring a locate to `want` against ring b?") takes one lane
//     per witness with an early exit by __any_sync; a single-witness locate splits the ring's edges over the lanes:
//     on-boundary = any(ballot), inside = popc(ballot(crossing)) & 1. Totals by warp sums, the count of rings with a
//     nonzero count by popc(ballot), and a `clean` flag (every total 0).
//  2. GridBox .. GridPairs: candidate edge pairs from a uniform grid per shape (g x g cells over the box of its clean
//     edges, g ~ sqrt(E)). Each edge goes into every cell its closed box covers; a pair of edges whose closed boxes
//     meet is classified only in the cell that holds the lower-left corner of the boxes' intersection (the reference
//     point), so it is counted once, by the same classify<T>(earlier edge, later edge) call as EdgePairs. Counts are
//     integer atomics: identical to EdgePairs for any order of the cells' entries.
//  3. BoundaryRows, EdtRowsH, EdtColsV: the 4-adjacent boundary (foreground with a 4-neighbour in the background, the
//     outside counts as background: objdel.reward.cache.boundary) by word shifts, and the exact squared Euclidean
//     distance transform, separable: row pass one warp per row (1D distances from clz / ctz on the packed words),
//     column pass one thread per column (Felzenszwalb-Huttenlocher's lower envelope of parabolas in Meijster's integer
//     form); both coalesced. The value min over sources of dx^2 + dy^2 is unique, so it is the same integer for every
//     correct algorithm.
//  4. Directed, Select: boundary F and HD95 of mask pairs from the EDT of each mask's boundary: per direction the number
//     of boundary pixels, the number within tau^2, and the order statistics at ranks k = floor(p (n - 1) / q) and k + 1
//     (bitwise radix select, one warp per direction, counts by warp sums).
//  5. SetInter, SetTerms: the intersection-count matrix of a set of masks (K samples then k raters; one warp per pair,
//     no atomics), then per set the squared GED (d = 1 - IoU) and HM-IoU by an exact integer dynamic program over the
//     rater quota states (each rater matched floor(K / k) or ceil(K / k) times), objective in 32.32 fixed point, ties
//     to the lexicographically smallest assignment.

// ===================================================================== warp helpers
using u64 = std::uint64_t;

struct Warp {
  int lane, size;  // GPU: lane in [0, 32), size 32. CPU build: lane 0, size 1
};

#ifdef __CUDACC__
template <class F>
__global__ void warp_kernel(i64 n, F f) {
  const i64 i = (i64(blockIdx.x) * blockDim.x + threadIdx.x) >> 5;  // warp-uniform: whole warps return together
  if (i < n) f(i, Warp{int(threadIdx.x & 31), 32});
}

template <class F>
void for_each_warp(Stream st, i64 n, const F& f) {
  constexpr int kThreads = 128;  // 4 warps per block
  if (n <= 0) return;
  warp_kernel<<<unsigned((n * 32 + kThreads - 1) / kThreads), kThreads, 0, st>>>(n, f);
  cuda_check(cudaGetLastError());
}
#else
template <class F>
void for_each_warp(Stream, i64 n, const F& f) {
  for (i64 i = 0; i < n; ++i) f(i, Warp{0, 1});
}
#endif

// Every lane of the warp must call these (warp-uniform control flow up to the call).
HD inline u32 wballot(bool b) {
#ifdef __CUDA_ARCH__
  return __ballot_sync(0xFFFFFFFFu, b);
#else
  return b ? 1u : 0u;
#endif
}
HD inline bool wany(bool b) { return wballot(b) != 0; }
HD inline bool wodd(bool b) { return (popc(wballot(b)) & 1) != 0; }  // parity of the lanes' flags
HD inline i64 wsum(i64 v) {  // integer: the same value for any order
#ifdef __CUDA_ARCH__
  for (int o = 16; o > 0; o >>= 1) v += __shfl_xor_sync(0xFFFFFFFFu, v, o);
#endif
  return v;
}
HD inline i64 wmin(i64 v) {
#ifdef __CUDA_ARCH__
  for (int o = 16; o > 0; o >>= 1) {
    const i64 u = __shfl_xor_sync(0xFFFFFFFFu, v, o);
    v = u < v ? u : v;
  }
#endif
  return v;
}
HD inline void wsync() {  // orders the lanes' global-memory writes before the next phase's reads
#ifdef __CUDA_ARCH__
  __syncwarp();
#endif
}

// ===================================================================== 1. validity totals, one warp per shape
struct ShapeTotalsWarp {
  Shapes s;
  Rings r;
  const u8* validate_wide;
  i64* defects;     // (n, 13), written whole
  u8* clean;        // per shape: every total is 0
  i32* bad_rings;   // per shape: rings with a nonzero count after the totals step

  // ShapeTotals::locate2 with the ring's edges split over the lanes.
  template <class T>
  HD int locate(i64 g, i64 qx, i64 qy, Warp w) const {
    const i32* box = r.box + 4 * g;
    if (qx < 2 * i64(box[0]) || qx > 2 * i64(box[1]) || qy < 2 * i64(box[2]) || qy > 2 * i64(box[3])) return -1;
    const i32* c = r.cxy + 2 * s.ring_start[g];
    const i32 n = r.len[g];
    bool on = false, in = false;
    for (i32 k = w.lane; k < n; k += w.size) {
      const i32 k1 = k + 1 == n ? 0 : k + 1;
      const i64 ax = 2 * i64(c[2 * k]), ay = 2 * i64(c[2 * k + 1]), bx = 2 * i64(c[2 * k1]), by = 2 * i64(c[2 * k1 + 1]);
      if (in_box(ax, ay, bx, by, qx, qy) && orient<T>(ax, ay, bx, by, qx, qy) == 0) on = true;
      else if ((ay <= qy && qy < by) || (by <= qy && qy < ay)) {
        const T num = T(ax - qx) * T(by - ay) + T(qy - ay) * T(bx - ax);
        if (by > ay ? num > 0 : num < 0) in = !in;
      }
    }
    if (wany(on)) return 0;  // ShapeTotals returns 0 at the first edge through the point, whatever came before
    return wodd(in) ? 1 : -1;
  }

  template <class T>
  HD int locate_region(i64 r0, i32 R, i32 o, i64 qx, i64 qy, Warp w) const {
    const int loc = locate<T>(r0 + o, qx, qy, w);
    if (loc == 0) return 0;
    bool in = loc > 0;
    for (i32 h = 0; h < R; ++h) {
      if (s.roles[r0 + h] != o || r.fewer[r0 + h]) continue;
      const int lh = locate<T>(r0 + h, qx, qy, w);
      if (lh == 0) return 0;
      if (lh > 0) in = !in;
    }
    return in ? 1 : -1;
  }

  // ShapeTotals::any_witness: the single witness by the edge-split locate, many witnesses one per lane (each lane runs
  // the serial locate), 32 at a time, stopping at the first chunk with a hit. "Any" does not depend on the order.
  template <class T>
  HD bool any_witness(const ShapeTotals& st, i64 r0, i32 R, i32 a, int want, bool single, i32 b, bool region,
                      Warp w) const {
    const i32* c = r.cxy + 2 * s.ring_start[r0 + a];
    const i32 n = r.len[r0 + a];
    if (single) {
      if (n <= 0) return false;
      const i64 x = 2 * i64(c[0]), y = 2 * i64(c[1]);
      return (region ? locate_region<T>(r0, R, b, x, y, w) : locate<T>(r0 + b, x, y, w)) == want;
    }
    for (i32 base = 0; base < 2 * n; base += w.size) {  // warp-uniform trip count
      const i32 i = base + w.lane;
      bool hit = false;
      if (i < 2 * n) {
        const i32 k = i >> 1, k1 = k + 1 == n ? 0 : k + 1;
        const i64 x = (i & 1) ? i64(c[2 * k]) + c[2 * k1] : 2 * i64(c[2 * k]);
        const i64 y = (i & 1) ? i64(c[2 * k + 1]) + c[2 * k1 + 1] : 2 * i64(c[2 * k + 1]);
        hit = (region ? st.locate2_region<T>(r0, R, b, x, y) : st.locate2<T>(r0 + b, x, y)) == want;
      }
      if (wany(hit)) return true;
    }
    return false;
  }

  template <class T>
  HD void run(i64 sh, Warp w) const {
    const ShapeTotals st{s, r, validate_wide, defects};
    const i64 r0 = s.shape_start[sh];
    const i32 R = i32(s.shape_start[sh + 1] - r0);
    const i32* roles = s.roles + r0;
    const i32 *pcross = r.pairs + sh * 2 * kPairs, *ptouch = pcross + kPairs;
    auto at = [&](i32 ring, int d) -> i32& { return r.counts[(r0 + ring) * kNumDefects + d]; };
    auto fewer = [&](i32 ring) { return r.fewer[r0 + ring] != 0; };
    auto apart = [&](i32 a, i32 b) { return pcross[a * kMaxRings + b] == 0 && ptouch[a * kMaxRings + b] == 0; };
    auto apart_region = [&](i32 a, i32 o) {
      if (!apart(a, o)) return false;
      for (i32 h = 0; h < R; ++h)
        if (roles[h] == o && !apart(a, h)) return false;
      return true;
    };
    // Phase 1, one lane per ring: ring-ring row sums, winding.
    for (i32 a = w.lane; a < R; a += w.size) {
      i32 c = 0, t = 0;
      for (i32 b = 0; b < R; ++b) c += pcross[a * kMaxRings + b], t += ptouch[a * kMaxRings + b];
      at(a, kRingRingCrossing) += c;
      at(a, kRingRingTouching) += t;
      if (!fewer(a)) {
        const i64 a2 = r.area2[r0 + a];
        if ((roles[a] < 0 && a2 < 0) || (roles[a] >= 0 && a2 > 0)) at(a, kWrongWinding) = 1;
      }
    }
    // Phase 2, the whole warp per test (uniform loops): containment by witness points; lane 0 writes.
    i64 nested_h = 0, nested_o = 0;
    for (i32 h = 0; h < R; ++h) {
      const i32 p = roles[h];
      if (p < 0 || fewer(h) || fewer(p)) continue;
      if (any_witness<T>(st, r0, R, h, -1, apart(h, p), p, false, w) && w.lane == 0) at(h, kHoleOutsideOuter) = 1;
    }
    for (i32 a = 0; a < R; ++a) {
      if (fewer(a)) continue;
      for (i32 b = a + 1; b < R; ++b) {
        if (fewer(b)) continue;
        if (roles[a] >= 0 && roles[a] == roles[b]) {
          const bool single = apart(a, b);
          if (any_witness<T>(st, r0, R, b, 1, single, a, false, w) ||
              any_witness<T>(st, r0, R, a, 1, single, b, false, w)) {
            if (w.lane == 0) at(a, kNestedHoles) += 1, at(b, kNestedHoles) += 1;
            ++nested_h;
          }
        } else if (roles[a] < 0 && roles[b] < 0) {
          if (any_witness<T>(st, r0, R, b, 1, apart_region(b, a), a, true, w) ||
              any_witness<T>(st, r0, R, a, 1, apart_region(a, b), b, true, w)) {
            if (w.lane == 0) at(a, kNestedOuters) += 1, at(b, kNestedOuters) += 1;
            ++nested_o;
          }
        }
      }
    }
    wsync();
    // Phase 3: totals by warp sums, bad rings by ballot + popc.
    i64* tot = defects + sh * kNumDefects;
    bool any = false;
    for (int d = 0; d < kNumDefects; ++d) {
      i64 v = 0;
      if (d == kRingRingCrossing || d == kRingRingTouching) {
        const i32* m = d == kRingRingCrossing ? pcross : ptouch;
        for (i32 q = w.lane; q < R * R; q += w.size)
          if (q / R < q % R) v += m[(q / R) * kMaxRings + q % R];
      } else if (d == kNestedHoles || d == kNestedOuters) {
        v = w.lane == 0 ? (d == kNestedHoles ? nested_h : nested_o) : 0;
      } else {
        for (i32 a = w.lane; a < R; a += w.size) v += at(a, d);
      }
      v = wsum(v);
      any = any || v != 0;
      if (w.lane == 0) tot[d] = v;
    }
    i32 bad = 0;
    for (i32 base = 0; base < R; base += w.size) {
      const i32 a = base + w.lane;
      bool nz = false;
      if (a < R)
        for (int d = 0; d < kNumDefects; ++d) nz = nz || at(a, d) != 0;
      bad += popc(wballot(nz));
    }
    if (w.lane == 0) clean[sh] = any ? 0 : 1, bad_rings[sh] = bad;
  }

  HD void operator()(i64 sh, Warp w) const {
    if (validate_wide[sh]) run<i128>(sh, w);
    else run<i64>(sh, w);
  }
};

// ===================================================================== 2. uniform-grid candidate edge pairs
// The clean edge starting at input vertex slot v (EdgePairs' item), if any.
HD inline bool slot_edge(const Shapes& s, const Rings& r, i64 v, i64& sh, i64& ring, i32& k) {
  const i64 nr = s.shape_start[s.n];
  ring = segment(s.ring_start, nr, v);
  k = i32(v - s.ring_start[ring]);
  if (r.fewer[ring] || k >= r.len[ring]) return false;
  sh = segment(s.shape_start, s.n, ring);
  return true;
}

struct GridView {
  const i32* g;           // cells per side per shape (host-chosen, 0 = no grid)
  const i64* cell_start;  // shape s = cells [cell_start[s], cell_start[s + 1]) = g^2 cells, row-major (cy g + cx)
  i32* gbox;              // x0, y0, x1, y1 over the shape's clean edges; x1 < x0 when it has none
  i32* cnt;               // entries per cell (GridCount), then the fill cursor (zeroed by GridScan)
  i64* off;               // first entry of each cell, ncells + 1, nondecreasing over all shapes
  i64* sh_entries;        // entries per shape (GridCount), zeroed
  const i64* sh_base;     // first entry of each shape (host prefix over the shapes that use the grid)
  const u8* use;          // the shape's pairs come from the grid (else EdgePairsRest)
  i64* entries;           // vertex slot per entry
  i32* ering;             // its ring (global index), so a pair test needs no search
  i64 ncells;
};

// Cell column of coordinate v in [lo, hi]: monotone, so an edge's cell range holds every cell of its closed box.
HD inline i64 grid_cell(i64 v, i64 lo, i64 hi, i32 g) { return (v - lo) * g / (hi - lo + 1); }

struct GridBox {
  Shapes s;
  Rings r;
  GridView gv;
  HD void operator()(i64 sh) const {
    i32 x0 = 1, y0 = 1, x1 = 0, y1 = 0;
    bool first = true;
    for (i64 q = s.shape_start[sh]; q < s.shape_start[sh + 1]; ++q) {
      if (r.fewer[q]) continue;
      const i32* b = r.box + 4 * q;  // x0, x1, y0, y1
      x0 = first ? b[0] : mn(x0, b[0]), x1 = first ? b[1] : mx(x1, b[1]);
      y0 = first ? b[2] : mn(y0, b[2]), y1 = first ? b[3] : mx(y1, b[3]);
      first = false;
    }
    i32* o = gv.gbox + 4 * sh;
    o[0] = x0, o[1] = y0, o[2] = x1, o[3] = y1;
  }
};

struct GridCount {
  Shapes s;
  Rings r;
  GridView gv;
  HD void operator()(i64 v) const {
    i64 sh = 0, ring = 0;
    i32 k = 0;
    if (!slot_edge(s, r, v, sh, ring, k) || gv.g[sh] == 0) return;
    const Edge e = make_edge(r, s, ring, i32(ring - s.shape_start[sh]), k);
    const i32 g = gv.g[sh];
    const i32* b = gv.gbox + 4 * sh;
    const i64 cx0 = grid_cell(e.xmin, b[0], b[2], g), cx1 = grid_cell(e.xmax, b[0], b[2], g);
    const i64 cy0 = grid_cell(e.ymin, b[1], b[3], g), cy1 = grid_cell(e.ymax, b[1], b[3], g);
    for (i64 cy = cy0; cy <= cy1; ++cy)
      for (i64 cx = cx0; cx <= cx1; ++cx) atomic_add(gv.cnt + gv.cell_start[sh] + cy * g + cx, 1);
    atomic_add(gv.sh_entries + sh, (cx1 - cx0 + 1) * (cy1 - cy0 + 1));
  }
};

struct GridScan {  // per shape: cell offsets from the shape's base, cursors zeroed
  Shapes s;
  GridView gv;
  HD void operator()(i64 sh) const {
    i64 run = gv.sh_base[sh];
    for (i64 c = gv.cell_start[sh]; c < gv.cell_start[sh + 1]; ++c) {
      gv.off[c] = run;
      if (gv.use[sh]) run += gv.cnt[c];
      gv.cnt[c] = 0;
    }
  }
};

struct GridFill {
  Shapes s;
  Rings r;
  GridView gv;
  HD void operator()(i64 v) const {
    i64 sh = 0, ring = 0;
    i32 k = 0;
    if (!slot_edge(s, r, v, sh, ring, k) || !gv.use[sh]) return;
    const Edge e = make_edge(r, s, ring, i32(ring - s.shape_start[sh]), k);
    const i32 g = gv.g[sh];
    const i32* b = gv.gbox + 4 * sh;
    const i64 cx0 = grid_cell(e.xmin, b[0], b[2], g), cx1 = grid_cell(e.xmax, b[0], b[2], g);
    const i64 cy0 = grid_cell(e.ymin, b[1], b[3], g), cy1 = grid_cell(e.ymax, b[1], b[3], g);
    for (i64 cy = cy0; cy <= cy1; ++cy)
      for (i64 cx = cx0; cx <= cx1; ++cx) {
        const i64 c = gv.cell_start[sh] + cy * g + cx, at = gv.off[c] + atomic_inc(gv.cnt + c);
        gv.entries[at] = v, gv.ering[at] = i32(ring);
      }
  }
};

// One item per entry: the entry's edge against every later entry of its cell. EdgePairs' counting, unchanged.
struct GridPairs {
  Shapes s;
  Rings r;
  GridView gv;
  const u8* validate_wide;
  template <class T>
  HD void count(i64 sh, i64 ring, const Edge& e, const Edge& f) const {
    const PairKind kind = classify<T>(e, f);
    if (kind == kNone) return;
    if (e.ring == f.ring) {
      const int d = kind == kCross ? kSelfCrossing : kind == kPinch ? kVertexPinch
                  : kind == kOverlap ? kCollinearOverlap : kSpike;
      atomic_add(r.counts + ring * kNumDefects + d, 1);
    } else {
      i32* m = r.pairs + sh * 2 * kPairs + (kind == kCross ? 0 : kPairs);
      atomic_add(m + e.ring * kMaxRings + f.ring, 1);
      atomic_add(m + f.ring * kMaxRings + e.ring, 1);
    }
  }
  HD void operator()(i64 p) const {
    const i64 c = segment(gv.off, gv.ncells, p), end = gv.off[c + 1];
    const i64 sh = segment(gv.cell_start, s.n, c), r0 = s.shape_start[sh];
    const i32 g = gv.g[sh];
    const i32* b = gv.gbox + 4 * sh;
    const i64 local = c - gv.cell_start[sh], ccx = local % g, ccy = local / g;
    const i64 v = gv.entries[p];
    const i64 rv = gv.ering[p];
    const Edge ev = make_edge(r, s, rv, i32(rv - r0), i32(v - s.ring_start[rv]));
    for (i64 q = p + 1; q < end; ++q) {
      const i64 u = gv.entries[q];
      const i64 ru = gv.ering[q];
      const Edge eu = make_edge(r, s, ru, i32(ru - r0), i32(u - s.ring_start[ru]));
      if (!(ev.xmin <= eu.xmax && eu.xmin <= ev.xmax && ev.ymin <= eu.ymax && eu.ymin <= ev.ymax)) continue;
      // the reference point: lower-left corner of the boxes' intersection; its cell owns the pair
      if (grid_cell(mx(ev.xmin, eu.xmin), b[0], b[2], g) != ccx || grid_cell(mx(ev.ymin, eu.ymin), b[1], b[3], g) != ccy)
        continue;
      const bool first = v < u;  // slots grow with (ring, k): EdgePairs' "edge, then a later edge"
      const Edge& e = first ? ev : eu;
      const Edge& f = first ? eu : ev;
      const i64 ring = first ? rv : ru;
      if (validate_wide[sh]) count<i128>(sh, ring, e, f);
      else count<i64>(sh, ring, e, f);
    }
  }
};

struct EdgePairsRest {  // EdgePairs for the shapes that do not use the grid
  Shapes s;
  Rings r;
  const u8* validate_wide;
  const u8* use;
  HD void operator()(i64 v) const {
    i64 sh = 0, ring = 0;
    i32 k = 0;
    if (!slot_edge(s, r, v, sh, ring, k) || use[sh]) return;
    EdgePairs{s, r, validate_wide}(v);
  }
};

// ===================================================================== 3. boundary and exact squared EDT
constexpr i32 kNoSource = 0x7FFFFFFF;  // squared distance when the mask has no source pixel

struct Frames {  // bit-packed masks (gpucheck layout) and their pixel arrays (i32, row-major, mask after mask)
  const u32* words;
  const i64 *word_start, *pix_start, *row_start, *col_start;
  const i32 *H, *W;
  i64 n;
};

struct BoundaryRows {  // item = row; out has the input's layout; nb[m] += boundary pixels (nb zeroed)
  Frames f;
  u32* out;
  i64* nb;
  HD void operator()(i64 i) const {
    const i64 m = segment(f.row_start, f.n, i);
    const i32 y = i32(i - f.row_start[m]), H = f.H[m], wpr = words_per_row(f.W[m]);
    const u32* row = f.words + f.word_start[m] + i64(y) * wpr;
    u32* o = out + f.word_start[m] + i64(y) * wpr;
    i64 cnt = 0;
    for (i32 j = 0; j < wpr; ++j) {
      const u32 cur = row[j], up = y > 0 ? row[j - wpr] : 0u, dn = y + 1 < H ? row[j + wpr] : 0u;
      const u32 left = (cur << 1) | (j > 0 ? row[j - 1] >> 31 : 0u);           // pixel x - 1 at bit x
      const u32 right = (cur >> 1) | (j + 1 < wpr ? row[j + 1] << 31 : 0u);   // pixel x + 1 (0 past W: precondition)
      const u32 b = cur & ~(up & dn & left & right);
      o[j] = b;
      cnt += popc(b);
    }
    if (cnt) atomic_add(nb + m, cnt);
  }
};

// Pass 1, one warp per row: g = distance along the row to the nearest source pixel, H + W when the row has none. Lane
// b takes pixel 32 j + b of word j (stores coalesced); the nearest source at or left of it is the highest set bit at
// or below b (clz), else the last source of the earlier words, carried warp-uniformly; the right side likewise (ctz).
struct EdtRowsH {
  Frames f;
  i32* g;
  HD void operator()(i64 i, Warp w) const {
    const i64 m = segment(f.row_start, f.n, i);
    const i32 y = i32(i - f.row_start[m]), H = f.H[m], W = f.W[m], wpr = words_per_row(W), inf = H + W;
    const u32* row = f.words + f.word_start[m] + i64(y) * wpr;
    i32* o = g + f.pix_start[m] + i64(y) * W;
    i64 last = -1;  // column of the last source in the words before j
    for (i32 j = 0; j < wpr; ++j) {
      const u32 bits = row[j];
      for (i32 b = w.lane; b < 32; b += w.size) {
        const i32 x = 32 * j + b;
        if (x >= W) continue;
        const u32 low = bits & ((2u << b) - 1u);  // bits 0..b (2u << 31 wraps to 0: all bits)
        const i64 src = low ? i64(32 * j + 31 - clz(low)) : last;
        o[x] = src >= 0 ? i32(x - src) : inf;
      }
      if (bits) last = 32 * j + 31 - clz(bits);
    }
    i64 next = -1;  // column of the first source in the words after j
    for (i32 j = wpr - 1; j >= 0; --j) {
      const u32 bits = row[j];
      for (i32 b = w.lane; b < 32; b += w.size) {  // the same lane wrote o[x] above: no sync needed
        const i32 x = 32 * j + b;
        if (x >= W) continue;
        const u32 high = bits & ~((1u << b) - 1u);  // bits b..31
        const i64 src = high ? i64(32 * j + ctz(high)) : next;
        if (src >= 0 && src - x < o[x]) o[x] = i32(src - x);
      }
      if (bits) next = 32 * j + ctz(bits);
    }
  }
};

// Pass 2, one thread per column: d(y) = min over rows i of (y - i)^2 + g(i)^2, the lower envelope of parabolas
// (Felzenszwalb and Huttenlocher 2012, Algorithm 1) in Meijster et al.'s integer form: Sep(i, u) = floor((u^2 - i^2 +
// g(u)^2 - g(i)^2) / (2 (u - i))) is the last row where parabola i is at least as low as parabola u. Adjacent threads
// take adjacent columns, so the loads of g at one row and the stores of d are coalesced. Scratch s (parabola rows) and
// t (first row of each) in the mask's pixel block, entry q of column x at q W + x. A column whose rows are all H + W
// belongs to a mask with no source pixel: kNoSource.
struct EdtColsV {
  Frames f;
  const i32* g;
  i32 *s, *t, *d;
  HD void operator()(i64 c) const {
    const i64 m = segment(f.col_start, f.n, c);
    const i32 x = i32(c - f.col_start[m]), H = f.H[m], W = f.W[m];
    const i64 base = f.pix_start[m] + x, inf = H + W;
    auto G = [&](i64 r) -> i64 { return g[base + r * W]; };
    auto S = [&](i64 q) -> i32& { return s[base + q * W]; };
    auto T = [&](i64 q) -> i32& { return t[base + q * W]; };
    bool any = false;
    for (i32 r = 0; r < H && !any; ++r) any = G(r) < inf;
    if (!any) {
      for (i32 r = 0; r < H; ++r) d[base + i64(r) * W] = kNoSource;
      return;
    }
    auto F = [&](i64 r, i64 i) { const i64 gi = G(i); return (r - i) * (r - i) + gi * gi; };
    auto sep = [&](i64 a, i64 u) {
      const i64 ga = G(a), gu = G(u);
      return floordiv<i64>(u * u - a * a + gu * gu - ga * ga, 2 * (u - a));
    };
    i32 q = 0;
    S(0) = 0, T(0) = 0;
    for (i32 u = 1; u < H; ++u) {
      while (q >= 0 && F(T(q), S(q)) > F(T(q), u)) --q;
      if (q < 0) {
        q = 0, S(0) = u;
      } else {
        const i64 v = 1 + sep(S(q), u);
        if (v < H) ++q, S(q) = u, T(q) = i32(v);
      }
    }
    for (i32 u = H - 1; u >= 0; --u) {
      d[base + i64(u) * W] = i32(F(u, S(q)));
      if (u == T(q)) --q;
    }
  }
};

// ===================================================================== 4. boundary F and HD95 per mask pair
struct Dirs {  // direction q = 2 pair + dir: points = boundary of src, distances to the boundary of dst
  const i64 *a, *b;          // per pair
  const i64* list_start;     // per direction (n_dir + 1): its distances in list
  const i64* drow_start;     // per direction (n_dir + 1): its rows (items of Directed)
  i64 n_dir;
  HD i64 src(i64 q) const { return (q & 1) ? b[q >> 1] : a[q >> 1]; }
  HD i64 dst(i64 q) const { return (q & 1) ? a[q >> 1] : b[q >> 1]; }
};

struct Directed {  // item = (direction, row): hits within tau^2 and the distances, appended in any order
  Frames f;
  Dirs dirs;
  const u32* bnd;
  const i32* d2;
  i64 tau2;
  i32* fill;   // per direction cursor (zeroed)
  i64* hit;    // per direction (zeroed)
  i32* list;
  HD void operator()(i64 i) const {
    const i64 q = segment(dirs.drow_start, dirs.n_dir, i), sm = dirs.src(q), dm = dirs.dst(q);
    const i32 y = i32(i - dirs.drow_start[q]), W = f.W[sm], wpr = words_per_row(W);
    const u32* row = bnd + f.word_start[sm] + i64(y) * wpr;
    const i32* dist = d2 + f.pix_start[dm] + i64(y) * W;
    i64 h = 0;
    for (i32 j = 0; j < wpr; ++j)
      for (u32 bits = row[j]; bits; bits &= bits - 1) {
        const i32 dv = dist[32 * j + ctz(bits)];
        h += dv <= tau2;
        list[dirs.list_start[q] + atomic_inc(fill + q)] = dv;
      }
    if (h) atomic_add(hit + q, h);
  }
};

// One warp per direction: the order statistics at ranks k = floor(p_num (n - 1) / p_den) and min(k + 1, n - 1). The
// k-th smallest value is the largest v with #(x < v) <= k, built bit by bit from the top (31 passes, each a warp sum).
struct Select {
  Dirs dirs;
  const i32* list;
  i64 p_num, p_den;
  i64 *rank, *lo, *hi;
  HD void operator()(i64 q, Warp w) const {
    const i32* x = list + dirs.list_start[q];
    const i64 n = dirs.list_start[q + 1] - dirs.list_start[q];
    if (n == 0) {
      if (w.lane == 0) rank[q] = -1, lo[q] = -1, hi[q] = -1;
      return;
    }
    const i64 k = p_num * (n - 1) / p_den;
    i64 ans = 0;
    for (int bit = 30; bit >= 0; --bit) {
      const i64 t = ans | (i64(1) << bit);
      i64 c = 0;
      for (i64 j = w.lane; j < n; j += w.size) c += x[j] < t;
      if (wsum(c) <= k) ans = t;
    }
    i64 le = 0, above = i64(1) << 40;
    for (i64 j = w.lane; j < n; j += w.size) {
      le += x[j] <= ans;
      if (x[j] > ans && x[j] < above) above = x[j];
    }
    le = wsum(le), above = wmin(above);
    if (w.lane == 0) rank[q] = k, lo[q] = ans, hi[q] = (k + 1 >= n || le >= k + 2) ? ans : above;
  }
};

// ===================================================================== 5. GED and HM-IoU over sets of masks
struct SetView {
  const u32* words;
  const i64* word_start;  // per mask
  const i32 *H, *W;       // per mask (one frame per set)
  const i64* mask_start;  // per set (n + 1)
  const i32* K;           // per set: samples first, then raters
  const i64* pair_start;  // per set: M (M + 1) / 2 pairs i <= j
  const i64* inter_start; // per set: M^2
  i64 n;
};

struct SetInter {  // one warp per pair (i <= j): |m_i & m_j| over the frame's words, no atomics
  SetView v;
  i64* inter;
  HD void operator()(i64 t, Warp w) const {
    const i64 s = segment(v.pair_start, v.n, t), m0 = v.mask_start[s];
    const i64 M = v.mask_start[s + 1] - m0;
    i64 local = t - v.pair_start[s], i = 0;
    while (local >= M - i) local -= M - i, ++i;
    const i64 j = i + local;
    const i64 nw = i64(v.H[m0]) * words_per_row(v.W[m0]);
    const u32 *a = v.words + v.word_start[m0 + i], *b = v.words + v.word_start[m0 + j];
    i64 c = 0;
    for (i64 k = w.lane; k < nw; k += w.size) c += popc(a[k] & b[k]);
    c = wsum(c);
    if (w.lane == 0) inter[v.inter_start[s] + i * M + j] = c, inter[v.inter_start[s] + j * M + i] = c;
  }
};

HD inline u64 iou_fixed(i64 inter, i64 a, i64 b) {  // floor(2^32 IoU), 2^32 for two empty masks
  const i64 uni = a + b - inter;
  return uni ? (u64(inter) << 32) / u64(uni) : (u64(1) << 32);
}

constexpr i64 kMaxHmStates = i64(1) << 16;

// Rater quota states of HM-IoU for K samples and k raters: (cap + 1)^k, cap = ceil(K / k); 0 when above the limit.
HD inline i64 hm_states(i64 K, i64 k) {
  const i64 cap = (K + k - 1) / k;
  i64 S = 1;
  for (i64 r = 0; r < k; ++r) {
    S *= cap + 1;
    if (S > kMaxHmStates) return 0;
  }
  return S;
}

struct SetTerms {  // per set, one thread: GED and HM-IoU from the integer matrix, in a fixed order
  SetView v;
  const i64* inter;
  const i64* state_start;  // per set: its HM states in best
  const i64* choice_start; // per set: K entries
  u64* best;
  i32* choice;
  u64* hm_q;
  double *hm_iou, *ged;
  HD void operator()(i64 s) const {
    const i64 M = v.mask_start[s + 1] - v.mask_start[s], K = v.K[s], k = M - K;
    const i64* X = inter + v.inter_start[s];
    auto I = [&](i64 i, i64 j) { return iou(X[i * M + j], X[i * M + i], X[j * M + j]); };
    double sxy = 0, sxx = 0, syy = 0;  // metrics._ged: 2 E d(S, Y) - E d(S, S') - E d(Y, Y'), d = 1 - IoU
    for (i64 i = 0; i < K; ++i)
      for (i64 j = K; j < M; ++j) sxy += 1.0 - I(i, j);
    for (i64 i = 0; i < K; ++i)
      for (i64 j = 0; j < K; ++j) sxx += 1.0 - I(i, j);
    for (i64 i = K; i < M; ++i)
      for (i64 j = K; j < M; ++j) syy += 1.0 - I(i, j);
    const double gd = 2.0 * (sxy / double(K * k)) - sxx / double(K * K) - syy / double(k * k);
    ged[s] = gd > 0.0 ? gd : 0.0;

    i32* ch = choice + choice_start[s];
    const i64 S = state_start[s + 1] - state_start[s];
    if (S == 0) {  // too many quota states for the device program
      for (i64 i = 0; i < K; ++i) ch[i] = -1;
      hm_q[s] = 0, hm_iou[s] = -1.0;
      return;
    }
    u64* B = best + state_start[s];
    const i64 base = K / k, cap = (K + k - 1) / k, R = cap + 1;
    const u64 kNeg = ~u64(0);
    auto Q = [&](i64 i, i64 r) { return iou_fixed(X[i * M + K + r], X[i * M + i], X[(K + r) * M + K + r]); };
    // B[state] = best total of the samples still to place (sum of counts .. K - 1); successors are larger states.
    for (i64 st = S - 1; st >= 0; --st) {
      i64 sum = 0, rest = st;
      bool full = true;
      for (i64 r = 0; r < k; ++r) {
        const i64 c = rest % R;
        rest /= R, sum += c, full = full && c >= base;
      }
      if (sum >= K) {
        B[st] = (sum == K && full) ? 0 : kNeg;
        continue;
      }
      u64 b = kNeg;
      i64 rad = 1;
      rest = st;
      for (i64 r = 0; r < k; ++r, rad *= R) {
        const i64 c = rest % R;
        rest /= R;
        if (c >= cap || B[st + rad] == kNeg) continue;
        const u64 val = Q(sum, r) + B[st + rad];
        if (b == kNeg || val > b) b = val;
      }
      B[st] = b;
    }
    i64 st = 0;
    double sum_iou = 0.0;
    for (i64 i = 0; i < K; ++i) {  // the smallest rater that keeps the optimum: lexicographically smallest assignment
      i64 rad = 1, rest = st;
      for (i64 r = 0; r < k; ++r, rad *= R) {
        const i64 c = rest % R;
        rest /= R;
        if (c >= cap || B[st + rad] == kNeg || Q(i, r) + B[st + rad] != B[st]) continue;
        ch[i] = i32(r);
        sum_iou += I(i, K + r);
        st += rad;
        break;
      }
    }
    hm_q[s] = B[0];
    hm_iou[s] = sum_iou / double(K);
  }
};
