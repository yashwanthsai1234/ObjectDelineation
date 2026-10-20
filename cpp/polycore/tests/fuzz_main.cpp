// Native fuzz / differential driver for polycore.hpp, built with ASan/UBSan by `make fuzz FUZZ_N=<shapes>`.
// Checks, on random and garbage input (int32-extreme coordinates, random grids, malformed offsets):
//  * no undefined behaviour or out-of-bounds access (ASan/UBSan abort otherwise);
//  * malformed or over-cap input is rejected by check_* with std::invalid_argument;
//  * rasterize == a brute-force per-pixel evaluation of the spec in __int128 (masks, per-ring winding
//    ranges, statistics), also on reused scratch;
//  * validate's sweep finds exactly the edge pairs of an all-pairs scan;
//  * score_batch results do not depend on the thread count.
#include <cstdio>
#include <cstdlib>
#include <random>
#include <stdexcept>

#include "polycore/polycore.hpp"

using namespace polycore;

#define CHECK(cond, ...)                         \
  do {                                           \
    if (!(cond)) {                               \
      std::printf("FAIL %s: ", #cond);           \
      std::printf(__VA_ARGS__);                  \
      std::printf("\n");                         \
      std::exit(1);                              \
    }                                            \
  } while (0)

// Winding of every centre from the spec, one ring or all rings (ring < 0), in __int128.
static std::vector<i32> brute_winding(const Rings& rs, i32 W, i32 H, i64 G, i32 only) {
  std::vector<i32> w(size_t(W) * H, 0);
  for (i32 ring = 0; ring < rs.nrings; ++ring) {
    if (only >= 0 && ring != only) continue;
    const i32 b = rs.begin(ring), n = rs.size(ring);
    for (i32 k = 0; k < n; ++k) {
      const i32 v0 = b + k, v1 = b + (k + 1 == n ? 0 : k + 1);
      const i128 x0 = i128(rs.x(v0)) * 2 * W, y0 = i128(rs.y(v0)) * 2 * H;
      const i128 x1 = i128(rs.x(v1)) * 2 * W, y1 = i128(rs.y(v1)) * 2 * H;
      for (i32 r = 0; r < H; ++r) {
        const i128 py = i128(2 * r + 1) * G;
        if (y0 == y1 || !(std::min(y0, y1) <= py && py < std::max(y0, y1))) continue;
        for (i32 c = 0; c < W; ++c) {
          const i128 num = (x0 - i128(2 * c + 1) * G) * (y1 - y0) + (py - y0) * (x1 - x0);
          if (y1 > y0 ? num > 0 : num < 0) w[size_t(r) * W + c] += y1 > y0 ? 1 : -1;
        }
      }
    }
  }
  return w;
}

template <class T>
static void brute_pairs(const ValidateScratch& s, ValidateResult& b) {
  const auto& E = s.edges;
  for (size_t i = 0; i < E.size(); ++i)
    for (size_t j = i + 1; j < E.size(); ++j)
      if (detail::boxes_meet(E[i], E[j])) detail::record_pair<T>(E[i], E[j], b);
}

template <class F>
static bool rejects(F f) {
  try {
    f();
  } catch (const std::invalid_argument&) {
    return true;
  }
  return false;
}

int main(int argc, char** argv) {
  const long iters = argc > 1 ? std::atol(argv[1]) : 20000;
  std::mt19937_64 rng(12345);
  auto U = [&](i64 lo, i64 hi) { return std::uniform_int_distribution<i64>(lo, hi)(rng); };
  const i32 picks[] = {INT32_MIN, INT32_MIN + 1, -1000000000, -1, 0, 1, 499, 500, 1000, 1000000000, INT32_MAX};
  RasterScratch rsc;
  ValidateScratch vsc;
  ValidateResult v;
  long checked = 0, rejected = 0;
  for (long it = 0; it < iters; ++it) {
    // ---- input: random offsets (sometimes malformed or over the caps) and coordinates
    const int mode = int(U(0, 3));
    const i32 R = i32(U(0, 12));
    std::vector<i32> off(R + 1, 0);
    for (i32 r = 0; r < R; ++r) off[r + 1] = off[r] + i32(U(0, U(0, 20) == 0 ? 1002 : 9));
    if (R && U(0, 20) == 0) off[U(0, R)] = i32(U(-3, 40));  // garbage offset
    const i32 nv = std::max<i32>(0, off[R] - (U(0, 20) == 0));  // sometimes one vertex short
    std::vector<i32> xy;
    for (i32 v_ = 0; v_ < 2 * nv; ++v_)
      xy.push_back(mode == 0 ? i32(U(0, 1000))
                 : mode == 1 ? i32(U(0, 20)) * 50
                 : mode == 2 ? picks[U(0, 10)]
                             : i32(U(INT32_MIN, INT32_MAX)));
    const Rings rs{xy.data(), off.data(), R};
    const i32 W = i32(U(1, 24)), H = i32(U(1, 24));
    const i64 G = U(0, 3) == 0 ? U(1, kMaxGrid) : 1000;
    bool ok = !rejects([&] { check_rings(rs, nv); });
    bool malformed = R > kMaxRings || off[0] < 0;
    for (i32 r = 0; r < R; ++r) malformed |= off[r + 1] < off[r] || off[r + 1] - off[r] > kMaxRingVertices;
    malformed |= off[R] > nv;
    CHECK(ok == !malformed, "check_rings accepted=%d malformed=%d it=%ld", ok, malformed, it);
    CHECK(rejects([&] { check_dims(0, H); }) && rejects([&] { check_dims(W, kMaxDim + 1); }), "dims");
    CHECK(rejects([&] { check_grid(0); }) && rejects([&] { check_grid(kMaxGrid + 1); }), "grid");
    if (!ok) {
      ++rejected;
      continue;
    }
    // ---- rasterize vs the spec (twice: the second call reuses the scratch)
    std::vector<u8> eo(size_t(W) * H, 7), nz(size_t(W) * H, 7);
    std::vector<i32> wmin(R), wmax(R);
    for (int rep = 0; rep < 2; ++rep) {
      const RasterStats st = rasterize(rs, W, H, G, rsc, RasterOut{eo.data(), nz.data(), wmin.data(), wmax.data()});
      const std::vector<i32> w = brute_winding(rs, W, H, G, -1);
      i64 aeo = 0, anz = 0, multi = 0;
      for (size_t p = 0; p < w.size(); ++p) {
        CHECK(eo[p] == u8(w[p] & 1) && nz[p] == u8(w[p] != 0), "mask it=%ld p=%zu", it, p);
        aeo += w[p] & 1, anz += w[p] != 0, multi += w[p] >= 2 || w[p] <= -2;
      }
      CHECK(st.area_evenodd == aeo && st.area_nonzero == anz && st.px_multi == multi, "stats it=%ld", it);
      i32 mixed = 0;
      for (i32 r = 0; r < R; ++r) {
        const std::vector<i32> wr = brute_winding(rs, W, H, G, r);
        const i32 lo = std::min(0, *std::min_element(wr.begin(), wr.end()));
        const i32 hi = std::max(0, *std::max_element(wr.begin(), wr.end()));
        CHECK(wmin[r] == lo && wmax[r] == hi, "ring range it=%ld ring=%d", it, r);
        mixed += lo < 0 && hi > 0;
      }
      CHECK(st.rings_mixed_sign == mixed, "mixed it=%ld", it);
    }
    // ---- validate: the sweep must see every pair an all-pairs scan sees
    std::vector<i32> roles(R, -1);
    for (i32 r = 1; r < R; ++r)
      if (U(0, 1)) roles[r] = 0;
    validate(rs, R ? roles.data() : nullptr, G, v, vsc);
    ValidateResult b;
    b.reset(R);
    if (detail::max_abs_coord(rs) < (i64(1) << 28)) brute_pairs<i64>(vsc, b);
    else brute_pairs<i128>(vsc, b);
    for (i32 r = 0; r < R; ++r)
      for (Defect d : {kSelfCrossing, kVertexPinch, kCollinearOverlap, kSpike})
        CHECK(v.at(r, d) == b.at(r, d), "sweep vs all pairs it=%ld ring=%d class=%d", it, r, int(d));
    CHECK(v.pair_cross == b.pair_cross && v.pair_touch == b.pair_touch, "pair matrices it=%ld", it);
    ++checked;
  }
  // ---- batch results independent of the thread count
  std::vector<std::vector<i32>> xys(64), offs(64);
  std::vector<Rollout> ros(64);
  for (int i = 0; i < 64; ++i) {
    offs[i] = {0, i32(U(3, 40))};
    for (i32 k = 0; k < 2 * offs[i][1]; ++k) xys[i].push_back(i32(U(0, 1000)));
    ros[i].rings = Rings{xys[i].data(), offs[i].data(), 1};
    ros[i].W = ros[i].H = i32(U(8, 128));
  }
  std::vector<RolloutScore> s1(64), s8(64);
  score_batch(ros.data(), 64, 1000, 1, s1.data());
  score_batch(ros.data(), 64, 1000, 8, s8.data());
  for (int i = 0; i < 64; ++i)
    CHECK(s1[i].defects == s8[i].defects && s1[i].raster.area_nonzero == s8[i].raster.area_nonzero &&
              s1[i].raster.area_evenodd == s8[i].raster.area_evenodd,
          "batch i=%d", i);
  std::printf("ok: %ld shapes checked against the spec and all-pairs scan, %ld malformed rejected, batch deterministic\n",
              checked, rejected);
  return 0;
}
