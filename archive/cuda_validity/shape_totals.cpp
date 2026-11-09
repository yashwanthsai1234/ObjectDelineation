// ShapeTotals of cuda/gpucheck.cpp at commit 1f64bd7 (D5), the last commit that ran it (archive/cuda_validity/README.md).
// Not built, not included.

// Ring-ring counts, winding and containment of one shape (polycore detail::validate_geometry, second half).
struct ShapeTotals {
  Shapes s;
  Rings r;
  const u8* validate_wide;
  i64* defects;

  // Point in doubled coordinates against clean ring g: -1 strictly outside, 0 on the boundary, 1 inside.
  template <class T>
  HD int locate2(i64 g, i64 qx, i64 qy) const {
    const i32* box = r.box + 4 * g;
    if (qx < 2 * i64(box[0]) || qx > 2 * i64(box[1]) || qy < 2 * i64(box[2]) || qy > 2 * i64(box[3])) return -1;
    const i32* c = r.cxy + 2 * s.ring_start[g];
    const i32 n = r.len[g];
    bool in = false;
    for (i32 k = 0; k < n; ++k) {
      const i32 k1 = k + 1 == n ? 0 : k + 1;
      const i64 ax = 2 * i64(c[2 * k]), ay = 2 * i64(c[2 * k + 1]), bx = 2 * i64(c[2 * k1]), by = 2 * i64(c[2 * k1 + 1]);
      if (in_box(ax, ay, bx, by, qx, qy) && orient<T>(ax, ay, bx, by, qx, qy) == 0) return 0;
      if ((ay <= qy && qy < by) || (by <= qy && qy < ay)) {
        const T num = T(ax - qx) * T(by - ay) + T(qy - ay) * T(bx - ax);
        if (by > ay ? num > 0 : num < 0) in = !in;
      }
    }
    return in ? 1 : -1;
  }

  // Region of outer ring o (local) with its holes, even-odd over all of them.
  template <class T>
  HD int locate2_region(i64 r0, i32 R, i32 o, i64 qx, i64 qy) const {
    const int loc = locate2<T>(r0 + o, qx, qy);
    if (loc == 0) return 0;
    bool in = loc > 0;
    for (i32 h = 0; h < R; ++h) {
      if (s.roles[r0 + h] != o || r.fewer[r0 + h]) continue;
      const int lh = locate2<T>(r0 + h, qx, qy);
      if (lh == 0) return 0;
      if (lh > 0) in = !in;
    }
    return in ? 1 : -1;
  }

  // Does a witness of ring a (vertex or edge midpoint) locate to `want` against ring b (or b's region)?
  // single: a neither crosses nor touches that boundary, so its first vertex decides.
  template <class T>
  HD bool any_witness(i64 r0, i32 R, i32 a, int want, bool single, i32 b, bool region) const {
    const i32* c = r.cxy + 2 * s.ring_start[r0 + a];
    const i32 n = r.len[r0 + a];
    auto loc = [&](i64 x, i64 y) { return region ? locate2_region<T>(r0, R, b, x, y) : locate2<T>(r0 + b, x, y); };
    if (single) return n > 0 && loc(2 * i64(c[0]), 2 * i64(c[1])) == want;
    for (i32 k = 0; k < n; ++k) {
      const i32 k1 = k + 1 == n ? 0 : k + 1;
      if (loc(2 * i64(c[2 * k]), 2 * i64(c[2 * k + 1])) == want) return true;
      if (loc(i64(c[2 * k]) + c[2 * k1], i64(c[2 * k + 1]) + c[2 * k1 + 1]) == want) return true;
    }
    return false;
  }

  template <class T>
  HD void run(i64 sh) const {
    const i64 r0 = s.shape_start[sh];
    const i32 R = i32(s.shape_start[sh + 1] - r0);
    const i32* roles = s.roles + r0;
    const i32 *pcross = r.pairs + sh * 2 * kPairs, *ptouch = pcross + kPairs;
    auto at = [&](i32 ring, int d) -> i32& { return r.counts[(r0 + ring) * kNumDefects + d]; };
    auto fewer = [&](i32 ring) { return r.fewer[r0 + ring] != 0; };
    auto apart = [&](i32 a, i32 b) { return pcross[a * kMaxRings + b] == 0 && ptouch[a * kMaxRings + b] == 0; };
    auto apart_region = [&](i32 a, i32 o) {  // a vs outer o and all of o's holes
      if (!apart(a, o)) return false;
      for (i32 h = 0; h < R; ++h)
        if (roles[h] == o && !apart(a, h)) return false;
      return true;
    };
    i64* tot = defects + sh * kNumDefects;
    for (i32 a = 0; a < R; ++a)
      for (i32 b = 0; b < R; ++b) {
        at(a, kRingRingCrossing) += pcross[a * kMaxRings + b];
        at(a, kRingRingTouching) += ptouch[a * kMaxRings + b];
      }
    for (i32 a = 0; a < R; ++a) {
      if (fewer(a)) continue;
      const i64 a2 = r.area2[r0 + a];
      if ((roles[a] < 0 && a2 < 0) || (roles[a] >= 0 && a2 > 0)) at(a, kWrongWinding) = 1;
    }
    for (i32 h = 0; h < R; ++h) {
      const i32 p = roles[h];
      if (p < 0 || fewer(h) || fewer(p)) continue;
      if (any_witness<T>(r0, R, h, -1, apart(h, p), p, false)) at(h, kHoleOutsideOuter) = 1;
    }
    for (i32 a = 0; a < R; ++a) {
      if (fewer(a)) continue;
      for (i32 b = a + 1; b < R; ++b) {
        if (fewer(b)) continue;
        if (roles[a] >= 0 && roles[a] == roles[b]) {  // sibling holes
          const bool single = apart(a, b);
          if (any_witness<T>(r0, R, b, 1, single, a, false) || any_witness<T>(r0, R, a, 1, single, b, false)) {
            at(a, kNestedHoles) += 1, at(b, kNestedHoles) += 1;
            tot[kNestedHoles] += 1;
          }
        } else if (roles[a] < 0 && roles[b] < 0) {  // two outers
          if (any_witness<T>(r0, R, b, 1, apart_region(b, a), a, true) ||
              any_witness<T>(r0, R, a, 1, apart_region(a, b), b, true)) {
            at(a, kNestedOuters) += 1, at(b, kNestedOuters) += 1;
            tot[kNestedOuters] += 1;
          }
        }
      }
    }
    for (i32 a = 0; a < R; ++a)
      for (int d = 0; d < kNumDefects; ++d)
        if (d != kRingRingCrossing && d != kRingRingTouching && d != kNestedHoles &&
            d != kNestedOuters)
          tot[d] += at(a, d);
    for (i32 a = 0; a < R; ++a)
      for (i32 b = a + 1; b < R; ++b) {
        tot[kRingRingCrossing] += pcross[a * kMaxRings + b];
        tot[kRingRingTouching] += ptouch[a * kMaxRings + b];
      }
  }

  HD void operator()(i64 sh) const {
    if (validate_wide[sh]) run<i128>(sh);
    else run<i64>(sh);
  }
};
