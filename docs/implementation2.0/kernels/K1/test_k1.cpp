// test_k1.cpp -- exact-equality tests of the K1 polygonizer.
// References:
//   (1) ref_polygonize below: a literal sequential port of objdel.checker.polygonize (edge list in boundary_edges
//       order, link by start-vertex lookup with the left turn at saddles, cycles, min-key start, ring order by min edge
//       index, corners) and of checker.bitquad (vertex codes). Independent of both K1 and gpucheck.
//   (2) gpucheck::check_masks (REPO/cuda/gpucheck.cpp, its CPU build), when compiled with -DWITH_GPUCHECK.
// Every field must be equal: quads, ring_start, pts_start, ring_area, pts. K1 runs with several seed_rows values
// and chunk sizes; all must give the same result.
//   usage: test_k1 [cases=3000] [seed=1]
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <stdexcept>

#include "k1poly.h"
#include "masks_gen.h"
#ifdef WITH_GPUCHECK
#include "gpucheck.h"
#endif

using namespace k1gen;
using k1poly::Result;

// ---------------------------------------------------------------- reference (literal port of polygonize.py)
static void ref_one(const Mask& mk, Result& out) {
  const i32 H = mk.H, W = mk.W;
  auto px = [&](i32 x, i32 y) -> int { return x >= 0 && y >= 0 && x < W && y < H && mk.px[size_t(y) * W + x]; };
  // quads (checker.bitquad)
  i64 q[5] = {0, 0, 0, 0, 0};
  for (i32 y = 0; y <= H; ++y)
    for (i32 x = 0; x <= W; ++x) {
      const int c = px(x - 1, y - 1) | px(x, y - 1) << 1 | px(x - 1, y) << 2 | px(x, y) << 3;
      const int b = __builtin_popcount(unsigned(c));
      if (c == 6 || c == 9) ++q[2];
      else if (b == 1) ++q[0];
      else if (b == 2) ++q[1];
      else if (b == 3) ++q[3];
      else if (b == 4) ++q[4];
    }
  out.quads.insert(out.quads.end(), q, q + 5);
  // boundary_edges: up (d 0), right (1), down (2), left (3), each in np.nonzero (row-major) order
  static const int DX[4] = {1, 0, -1, 0}, DY[4] = {0, 1, 0, -1};
  std::vector<i32> ex, ey, ed;
  const int sx[4] = {0, 1, 0, -1}, sy[4] = {-1, 0, 1, 0};  // neighbour that must be background
  const int ox[4] = {0, 1, 1, 0}, oy[4] = {0, 0, 1, 1};   // start-vertex offset
  for (int d = 0; d < 4; ++d)
    for (i32 y = 0; y < H; ++y)
      for (i32 x = 0; x < W; ++x)
        if (px(x, y) && !px(x + sx[d], y + sy[d])) ex.push_back(x + ox[d]), ey.push_back(y + oy[d]), ed.push_back(d);
  const i64 E = i64(ex.size());
  const i64 base = out.ring_area.size();
  if (E == 0) {
    out.ring_start.push_back(base);
    return;
  }
  // link: edges by start vertex (stable, i.e. by index)
  std::map<i64, std::vector<i64>> by_start;
  for (i64 e = 0; e < E; ++e) by_start[i64(ey[e]) * (W + 1) + ex[e]].push_back(e);
  std::vector<i64> succ(E);
  for (i64 e = 0; e < E; ++e) {
    const i64 key = i64(ey[e] + DY[ed[e]]) * (W + 1) + (ex[e] + DX[ed[e]]);
    const auto& c = by_start.at(key);
    if (c.size() == 1) succ[e] = c[0];
    else if (c.size() == 2) {
      const int want = (ed[e] + 3) % 4;
      if (ed[c[0]] == want) succ[e] = c[0];
      else if (ed[c[1]] == want) succ[e] = c[1];
      else throw std::logic_error("ref: no left turn at saddle");
    } else throw std::logic_error("ref: vertex with 0 or > 2 outgoing edges");
  }
  // rings in order of their minimum edge index = order of first discovery in index order
  std::vector<char> seen(E, 0);
  for (i64 e0 = 0; e0 < E; ++e0) {
    if (seen[e0]) continue;
    std::vector<i64> cyc;
    for (i64 e = e0; !seen[e]; e = succ[e]) seen[e] = 1, cyc.push_back(e);
    const i64 L = i64(cyc.size());
    i64 best = 0, area2 = 0;
    auto key = [&](i64 e) { return (i64(ey[e]) * (W + 1) + ex[e]) * 4 + ed[e]; };
    for (i64 k = 0; k < L; ++k) {
      if (key(cyc[k]) < key(cyc[best])) best = k;
      const i64 e = cyc[k], x1 = ex[e] + DX[ed[e]], y1 = ey[e] + DY[ed[e]];
      area2 += i64(ex[e]) * y1 - x1 * ey[e];
    }
    for (i64 k = 0; k < L; ++k) {
      const i64 e = cyc[(best + k) % L], p = cyc[(best + k + L - 1) % L];
      if (ed[e] != ed[p]) out.pts.push_back(ex[e]), out.pts.push_back(ey[e]);
    }
    out.pts_start.push_back(i64(out.pts.size() / 2));
    out.ring_area.push_back(area2 / 2);
  }
  out.ring_start.push_back(i64(out.ring_area.size()));
}

static Result ref_polygonize(const std::vector<const Mask*>& ms) {
  Result out;
  out.ring_start = {0};
  out.pts_start = {0};
  for (const Mask* m : ms) ref_one(*m, out);
  return out;
}

// ---------------------------------------------------------------- comparison
template <class A, class B>
static bool same_vec(const A& a, const B& b, const char* what, std::string& why) {
  if (a.size() != b.size()) {
    why = std::string(what) + ": size " + std::to_string(a.size()) + " vs " + std::to_string(b.size());
    return false;
  }
  for (size_t k = 0; k < a.size(); ++k)
    if (i64(a[k]) != i64(b[k])) {
      why = std::string(what) + "[" + std::to_string(k) + "]: " + std::to_string(i64(a[k])) + " vs " +
            std::to_string(i64(b[k]));
      return false;
    }
  return true;
}

template <class R1, class R2>
static bool same(const R1& a, const R2& b, std::string& why) {
  return same_vec(a.quads, b.quads, "quads", why) && same_vec(a.ring_start, b.ring_start, "ring_start", why) &&
         same_vec(a.pts_start, b.pts_start, "pts_start", why) && same_vec(a.ring_area, b.ring_area, "ring_area", why) &&
         same_vec(a.pts, b.pts, "pts", why);
}

static void dump(const Mask& m) {
  std::fprintf(stderr, "mask %s %dx%d\n", m.name.c_str(), m.H, m.W);
  if (m.H * m.W > 4096) return;
  for (i32 y = 0; y < m.H; ++y) {
    for (i32 x = 0; x < m.W; ++x) std::fputc(m.px[size_t(y) * m.W + x] ? '#' : '.', stderr);
    std::fputc('\n', stderr);
  }
}

static int failures = 0;
static i64 checked = 0, rings_seen = 0;

static void check_batch(const std::vector<const Mask*>& ms, const char* label) {
  Batch b;
  for (const Mask* m : ms) b.add(*m);
  const auto in = b.view();
  const Result ref = ref_polygonize(ms);
  std::string why;
#ifdef WITH_GPUCHECK
  const gpucheck::Masks gin{in.n, in.bits, in.start, in.H, in.W};
  const gpucheck::MaskResult gref = gpucheck::check_masks(gin);
  if (!same(gref, ref, why)) {
    std::fprintf(stderr, "[%s] REFERENCES DISAGREE (gpucheck vs ref): %s\n", label, why.c_str());
    if (ms.size() == 1) dump(*ms[0]);
    ++failures;
    return;
  }
#endif
  const int seeds[] = {0, 1, 2, 3, 16};
  for (int K : seeds)
    for (int cfg = 0; cfg < 3; ++cfg) {  // one chunk; chunks of <= 700 pixels; scratch limit 1 byte (halving to one mask)
      const i64 chunk = cfg == 1 ? 700 : i64(1) << 25;
      k1poly::Options opt;
      opt.seed_rows = K, opt.chunk_pixels = chunk, opt.max_scratch_bytes = cfg == 2 ? 1 : i64(1) << 30;
      Result got;
      try {
        got = k1poly::polygonize_masks(in, opt);
      } catch (const std::exception& e) {
        why = std::string("exception: ") + e.what();
        got = Result{};
      }
      if (!why.empty() || !same(got, ref, why)) {
        std::fprintf(stderr, "[%s] MISMATCH seed_rows=%d cfg=%d: %s\n", label, K, cfg, why.c_str());
        if (ms.size() == 1) dump(*ms[0]);
        ++failures;
        return;
      }
    }
  checked += i64(ms.size());
  rings_seen += i64(ref.ring_area.size());
}

int main(int argc, char** argv) {
  const int cases = argc > 1 ? std::atoi(argv[1]) : 3000;
  const unsigned long long seed = argc > 2 ? std::strtoull(argv[2], nullptr, 10) : 1ull;
  std::mt19937_64 rng(seed);
  const auto t0 = std::chrono::steady_clock::now();
  std::printf("k1poly test: %s build, %d random cases, seed %llu%s\n", k1poly::cuda_build() ? "CUDA" : "CPU", cases,
              seed,
#ifdef WITH_GPUCHECK
              ", references: ref_polygonize + gpucheck::check_masks (CPU)"
#else
              ", reference: ref_polygonize"
#endif
  );

  // 1. adversarial, single masks
  std::vector<Mask> adv;
  for (auto [h, w] : std::vector<std::pair<int, int>>{{1, 1}, {1, 7}, {7, 1}, {2, 2}, {5, 33}, {33, 5}, {31, 32},
                                                      {32, 33}, {64, 64}, {65, 63}, {40, 97}}) {
    adv.push_back(full(h, w, 0));
    adv.push_back(full(h, w, 1));
    adv.push_back(checker(h, w));
    adv.push_back(concentric(h, w));
    adv.push_back(comb(h, w));
    adv.push_back(spiral(h, w));
    adv.push_back(diagonals(h, w));
    adv.push_back(stripes(h, w, true));
    adv.push_back(stripes(h, w, false));
    adv.push_back(frame_with_holes(h, w));
  }
  adv.push_back(blank(0, 5, "H=0"));
  adv.push_back(blank(5, 0, "W=0"));
  adv.push_back(blank(0, 0, "0x0"));
  adv.push_back(checker(256, 256));
  adv.push_back(comb(200, 256));
  adv.push_back(spiral(256, 256));
  adv.push_back(concentric(256, 256));
  adv.push_back(frame_with_holes(128, 160));
  for (int k = 0; k < 6; ++k) adv.push_back(blobs(rng, 512, 512, 3 + 4 * k, 20.0 + 15 * k));
  for (const Mask& m : adv) check_batch({&m}, m.name.c_str());
  std::printf("adversarial: %zu masks done, failures so far %d\n", adv.size(), failures);

  // 2. random single masks: sizes 1..80, densities 0.05..0.95; and blobs
  std::uniform_int_distribution<int> dim(1, 80);
  std::uniform_real_distribution<double> dens(0.05, 0.95);
  for (int c = 0; c < cases && failures < 5; ++c) {
    Mask m = (c % 4 == 3) ? blobs(rng, dim(rng) + 20, dim(rng) + 20, 1 + c % 7, 3.0 + c % 13) : noise(rng, dim(rng), dim(rng), dens(rng));
    check_batch({&m}, m.name.c_str());
  }
  std::printf("random singles: done, failures so far %d\n", failures);

  // 3. random batches of mixed sizes (chunking, empty masks inside a batch, per-mask ring offsets)
  for (int c = 0; c < std::max(1, cases / 20) && failures < 5; ++c) {
    std::vector<Mask> ms;
    const int nm = 1 + c % 17;
    for (int k = 0; k < nm; ++k) {
      const int kind = int(rng() % 5);
      const int h = int(rng() % 50), w = int(rng() % 50);
      if (kind == 0) ms.push_back(blank(h, w, "blank"));
      else if (kind == 1) ms.push_back(checker(h, w));
      else ms.push_back(noise(rng, h, w, dens(rng)));
    }
    std::vector<const Mask*> p;
    for (const Mask& m : ms) p.push_back(&m);
    check_batch(p, "batch");
  }
  const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::printf("%s: %lld masks, %lld rings checked, %d failures, %.1f s\n", failures ? "FAIL" : "PASS",
              (long long)checked, (long long)rings_seen, failures, s);
  return failures ? 1 : 0;
}
