// K3 CPU tests: the shared formulas and the CPU reference, with no GPU. Build and run: make cpu && ./build/test_cpu
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

#include "k3_ref.h"

using namespace k3;
using namespace k3ref;

static int g_fail = 0;
#define CHECK(cond, ...)                         \
  do {                                           \
    if (!(cond)) {                               \
      ++g_fail;                                  \
      std::printf("FAIL %s:%d: ", __FILE__, __LINE__); \
      std::printf(__VA_ARGS__);                  \
      std::printf("\n");                         \
    }                                            \
  } while (0)

static double ulp_err(float got, double want) {
  if (std::isnan(want)) return std::isnan(got) ? 0 : 1e9;
  if (std::isinf(want)) return (std::isinf(got) && (got > 0) == (want > 0)) ? 0 : 1e9;
  const float wf = (float)want;
  if (std::isinf(wf)) return (std::isinf(got) && (got > 0) == (wf > 0)) ? 0 : 1e9;  // rounds to +-inf
  float ulp = std::nextafter(std::fabs(wf), INFINITY) - std::fabs(wf);
  if (wf == 0.0f || !std::isfinite(ulp)) ulp = std::nextafter(0.0f, 1.0f);
  return std::fabs((double)got - want) / ulp;
}

static void test_exp_log() {
  std::mt19937 rng(1);
  double emax = 0, lmax = 0;
  float ex_at = 0, lx_at = 0;
  // dense sweep of exp's domain + random bit patterns for log
  for (int i = 0; i <= 2000000; ++i) {
    float x = -103.0f + 192.0f * (float)i / 2000000.0f;
    double e = ulp_err(k3_expf(x), std::exp((double)x));
    if (e > emax) { emax = e; ex_at = x; }
  }
  std::uniform_int_distribution<uint32_t> bits(0x00000001u, 0x7f7fffffu);
  for (int i = 0; i < 2000000; ++i) {
    float x = u2f(bits(rng));
    double e = ulp_err(k3_logf(x), std::log((double)x));
    if (e > lmax) { lmax = e; lx_at = x; }
  }
  std::printf("k3_expf max error %.3f ulp (at %.9g); k3_logf max error %.3f ulp (at %.9g)\n", emax, ex_at, lmax, lx_at);
  CHECK(emax <= 2.0, "exp error %.3f ulp", emax);
  CHECK(lmax <= 2.0, "log error %.3f ulp", lmax);
  CHECK(k3_expf(0.0f) == 1.0f, "exp(0)");
  CHECK(k3_expf(-K3_INF) == 0.0f, "exp(-inf)");
  CHECK(k3_expf(K3_INF) == K3_INF, "exp(inf)");
  CHECK(k3_expf(89.0f) == K3_INF, "exp overflow");
  CHECK(std::isnan(k3_expf(NAN)), "exp(nan)");
  CHECK(k3_logf(1.0f) == 0.0f, "log(1)");
  CHECK(k3_logf(0.0f) == -K3_INF, "log(0)");
  CHECK(std::isnan(k3_logf(-1.0f)), "log(-1)");
  CHECK(ulp_err(k3_logf(1e-45f), std::log((double)1e-45f)) <= 2, "log(subnormal)");
  CHECK(ulp_err(k3_logf(3e-39f), std::log((double)3e-39f)) <= 2, "log(subnormal 2)");
  CHECK(ulp_err(k3_expf(-100.0f), std::exp(-100.0)) <= 2, "exp -> subnormal");
  // bf16 rounding: ties to even and the carry into the exponent
  CHECK(f32_to_bf16_rne(u2f(0x3f808000u)) == 0x3f80, "tie to even (down)");
  CHECK(f32_to_bf16_rne(u2f(0x3f818000u)) == 0x3f82, "tie to even (up)");
  CHECK(f32_to_bf16_rne(u2f(0x3f80ffffu)) == 0x3f81, "round up");
  CHECK(f32_to_bf16_rne(u2f(0x7f7fffffu)) == 0x7f80, "overflow to inf");
}

static void test_preprocess() {
  std::mt19937 rng(2);
  PreParams p;  // Qwen3-VL: 512 x 512, ps 16, ms 2, tps 2, mean' = std' = 127.5
  const PatchGeom g = geom(p);
  CHECK(g.G == 1024 && g.F == 1536, "geometry %d x %d", g.G, g.F);
  const size_t img = (size_t)p.H * p.W * 3;
  const int n_src = 5;
  std::vector<uint8_t> src(img * n_src);
  for (auto& v : src) v = (uint8_t)(rng() & 0xff);
  // adversarial images: 0, 255, and every value 0..255 in every channel at varied positions
  std::fill(src.begin(), src.begin() + img, 0);
  std::fill(src.begin() + img, src.begin() + 2 * img, 255);
  for (size_t i = 0; i < img; ++i) src[2 * img + i] = (uint8_t)((i * 7 + i / 3) & 0xff);
  const std::vector<int> idx = {4, 0, 0, 1, 2, 3, 2, 4};  // duplicated rows, as in a generation batch
  const int n_out = (int)idx.size();
  std::vector<float> a((size_t)n_out * g.G * g.F), b(a.size());
  for (int hwc = 1; hwc >= 0; --hwc) {
    p.hwc = hwc;
    preprocess_ref(src.data(), idx.data(), n_out, p, a.data());
    preprocess_emul(src.data(), idx.data(), n_out, p, b.data());
    CHECK(std::memcmp(a.data(), b.data(), a.size() * 4) == 0, "preprocess emul != ref (hwc=%d)", hwc);
  }
  // all 256 values per channel: the formula against an independent evaluation in double, rounded once
  int bad = 0;
  for (int v = 0; v < 256; ++v) {
    float y = k3_norm((uint8_t)v, 127.5f, 127.5f);
    double want = ((double)v - 127.5) / 127.5;  // sub is exact here, so one rounding of the exact quotient
    if (y != (float)want) ++bad;
  }
  CHECK(bad == 0, "%d of 256 values differ from the correctly rounded quotient", bad);
  // a non-default geometry (non-square, 3 x 2 merged blocks) and other mean/std
  PreParams q;
  q.H = 96; q.W = 64; q.mean[0] = 122.77f; q.mean[1] = 116.746f; q.mean[2] = 104.094f;
  q.std_[0] = 68.5f; q.std_[1] = 66.63f; q.std_[2] = 70.32f;
  const PatchGeom gq = geom(q);
  std::vector<uint8_t> s2((size_t)q.H * q.W * 3 * 2);
  for (auto& v : s2) v = (uint8_t)(rng() & 0xff);
  std::vector<float> c((size_t)2 * gq.G * gq.F), d(c.size());
  preprocess_ref(s2.data(), nullptr, 2, q, c.data());
  preprocess_emul(s2.data(), nullptr, 2, q, d.data());
  CHECK(std::memcmp(c.data(), d.data(), c.size() * 4) == 0, "preprocess emul != ref (96 x 64)");
  std::printf("preprocess: emul == ref on %d + 2 images, both layouts\n", n_out);
}

static std::vector<int64_t> synth_prompt(int pre, int lh, int lw, int post, std::vector<int64_t>& types) {
  std::vector<int64_t> ids;
  types.clear();
  for (int i = 0; i < pre; ++i) { ids.push_back(1000 + i); types.push_back(0); }
  ids.push_back(151652); types.push_back(0);  // vision_start
  for (int i = 0; i < lh * lw; ++i) { ids.push_back(151655); types.push_back(1); }
  ids.push_back(151653); types.push_back(0);  // vision_end
  for (int i = 0; i < post; ++i) { ids.push_back(2000 + i); types.push_back(0); }
  return ids;
}

static void test_rope_case(const char* name, int pre, Grid grid, int post, int B, int W, std::mt19937& rng) {
  const int merge = 2;
  std::vector<int64_t> ptypes;
  auto pids = synth_prompt(pre, (int)(grid[1] / merge), (int)(grid[2] / merge), post, ptypes);
  const int Lp = (int)pids.size();
  RopeTemplate t = rope_template(pids.data(), ptypes.data(), Lp, {grid}, merge);
  std::vector<int64_t> ids((size_t)B * W), mask(ids.size()), types(ids.size());
  const int P = Lp + 7;  // prompt width with up to 7 left pads
  for (int b = 0; b < B; ++b) {
    int lpad = (int)(rng() % 8), L = (int)(rng() % (W - P + 1));
    if (b == 0) L = 0;              // empty completion
    if (b == 1) L = W - P;          // full width
    for (int j = 0; j < W; ++j) {
      size_t o = (size_t)b * W + j;
      ids[o] = 0; mask[o] = 0; types[o] = 0;
      if (j >= lpad && j < lpad + Lp) { ids[o] = pids[j - lpad]; types[o] = ptypes[j - lpad]; mask[o] = 1; }
      if (j >= P && j < P + L) { ids[o] = 3 + (int64_t)(rng() % 150000); mask[o] = 1; }
      if (j < P && j >= lpad + Lp) { ids[o] = 151643; mask[o] = 0; }  // pads between prompt and answer
    }
    if (b == 2 && L > 10) mask[(size_t)b * W + P + 5] = 0;  // a hole inside the answer
  }
  std::vector<int64_t> pe((size_t)3 * B * W), pr((size_t)3 * W);
  int err = -1;
  rope_emul(B, W, ids.data(), mask.data(), types.data(), t, pe.data(), &err);
  CHECK(err == 0, "%s: err %d", name, err);
  int bad = 0;
  for (int b = 0; b < B; ++b) {
    rope_ref_row(ids.data() + (size_t)b * W, mask.data() + (size_t)b * W, types.data() + (size_t)b * W, W, {grid},
                 merge, pr.data());
    for (int a = 0; a < 3; ++a)
      for (int j = 0; j < W; ++j) bad += pe[(size_t)a * B * W + (size_t)b * W + j] != pr[(size_t)a * W + j];
  }
  CHECK(bad == 0, "%s: %d positions differ", name, bad);
  // adversarial: wrong prompt token, image token in the answer, all-masked row
  auto ids2 = ids; ids2[(size_t)3 * W + 7 + 2] = 77;  // row 3, inside the prompt (lpad <= 7 < 7 + 2 ... in prompt)
  rope_emul(B, W, ids2.data(), mask.data(), types.data(), t, pe.data(), &err);
  CHECK(err & 1, "%s: prompt mismatch not flagged (err %d)", name, err);
  auto ty2 = types; ty2[(size_t)1 * W + W - 1] = 1;  // row 1 has a full-width answer
  rope_emul(B, W, ids.data(), mask.data(), ty2.data(), t, pe.data(), &err);
  CHECK(err & 2, "%s: image token in answer not flagged (err %d)", name, err);
  auto m2 = mask; for (int j = 0; j < W; ++j) m2[(size_t)4 * W + j] = 0;
  rope_emul(B, W, ids.data(), m2.data(), types.data(), t, pe.data(), &err);
  CHECK(err & 4, "%s: short row not flagged (err %d)", name, err);
  std::printf("rope %-22s Lp %d end %lld: emul == HF-port ref on %d rows x %d\n", name, Lp, (long long)t.end, B, W);
}

static void test_rope() {
  std::mt19937 rng(3);
  test_rope_case("qwen3vl 512x512", 3, {1, 32, 32}, 60, 16, 1900, rng);
  test_rope_case("non-square 1x32x64", 5, {1, 32, 64}, 9, 8, 900, rng);
}

static void fill_bf16(std::vector<uint16_t>& v, std::mt19937& rng, float scale) {
  std::normal_distribution<float> nd(0.0f, scale);
  for (auto& x : v) x = f32_to_bf16_rne(nd(rng));
}

static void test_lse() {
  std::mt19937 rng(4);
  const int N = 70, V = 1000, D = 64;
  std::vector<uint16_t> H((size_t)N * D), W((size_t)V * D);
  fill_bf16(H, rng, 1.0f);
  fill_bf16(W, rng, 0.4f);
  for (int k = 0; k < D; ++k) {
    W[(size_t)17 * D + k] = f32_to_bf16_rne(bf16_to_f32(H[(size_t)5 * D + k]) * 3.0f);  // row 5: peaked on 17
    H[(size_t)6 * D + k] = 0;                                                        // row 6: flat (all z = 0)
  }
  std::vector<int32_t> tgt(N);
  for (int n = 0; n < N; ++n) tgt[n] = (int32_t)(rng() % V);
  tgt[0] = 0; tgt[1] = V - 1; tgt[2] = -100; tgt[3] = V; tgt[5] = 17;
  std::vector<double> l64(N), e64(N);
  lse_ref_f64(H.data(), W.data(), tgt.data(), N, V, D, 1.0, l64.data(), e64.data());
  for (int S : {1, 3, 8, 16}) {
    std::vector<float> l(N), e(N);
    lse_ref(H.data(), W.data(), tgt.data(), N, V, D, 1.0f, S, l.data(), e.data());
    double dl = 0, de = 0;
    for (int n = 0; n < N; ++n) {
      dl = std::max(dl, std::fabs(l[n] - l64[n]));
      de = std::max(de, std::fabs(e[n] - e64[n]));
    }
    std::printf("lse_ref S=%-2d: max |logp - fp64| %.2e, max |entropy - fp64| %.2e\n", S, dl, de);
    CHECK(dl < 2e-5 && de < 2e-5, "lse_ref S=%d far from fp64 (%.2e, %.2e)", S, dl, de);
    CHECK(l[2] == 0.0f && l[3] == 0.0f, "invalid targets must give 0");
    CHECK(std::fabs(l[6] + std::log((double)V)) < 1e-5, "flat row logp %.7f", l[6]);
    CHECK(l[5] > -1e-3f, "peaked row logp %.7f", l[5]);
  }
  // temperature
  std::vector<float> lt(N);
  std::vector<double> lt64(N);
  lse_ref(H.data(), W.data(), tgt.data(), N, V, D, 1.0f / 0.7f, 4, lt.data(), nullptr);
  lse_ref_f64(H.data(), W.data(), tgt.data(), N, V, D, (double)(1.0f / 0.7f), lt64.data(), nullptr);
  double dt = 0;
  for (int n = 0; n < N; ++n) dt = std::max(dt, std::fabs(lt[n] - lt64[n]));
  CHECK(dt < 2e-5, "temperature 0.7: %.2e", dt);
}

static void test_is() {
  std::mt19937 rng(5);
  std::normal_distribution<float> nd(0.0f, 0.05f);
  const int B = 37, T = 700;
  std::vector<float> lp((size_t)B * T), sp(lp.size());
  std::vector<uint8_t> m(lp.size());
  for (size_t i = 0; i < lp.size(); ++i) {
    lp[i] = -std::fabs(nd(rng)) * 20;
    sp[i] = lp[i] + nd(rng);
    m[i] = (rng() % 10) != 0;
  }
  sp[3] = NAN;              // vLLM could not score
  m[4] = 0; sp[4] = NAN;    // masked NaN: (x * 0) = NaN -> 0
  sp[5] = -INFINITY;        // diff +inf -> FLT_MAX -> exp = inf
  lp[6] = -80.0f; sp[6] = 0; // diff -80 -> ratio ~ 0
  for (int j = 0; j < T; ++j) m[(size_t)7 * T + j] = 0;  // row 7 fully masked
  const float INF = K3_INF;
  struct Case { int mode; float lo, hi; } cases[] = {
      {TOKEN_TRUNCATE, -INF, 3.0f}, {TOKEN_MASK, 0.0f, 3.0f}, {TOKEN_TRUNCATE, 0.5f, 2.0f},
      {SEQUENCE_TRUNCATE, -INF, 3.0f}, {SEQUENCE_MASK, -INF, 3.0f}, {TOKEN_MASK, -INF, INF}};
  for (const Case& c : cases) {
    const bool seq = c.mode >= SEQUENCE_TRUNCATE;
    std::vector<float> r(seq ? B : (size_t)B * T), rs((size_t)B * IS_NSTATS), st(IS_NSTATS);
    is_ref(lp.data(), sp.data(), m.data(), B, T, c.mode, c.lo, c.hi, r.data(), rs.data(), st.data());
    // independent check in double with torch semantics
    int bad = 0;
    double dsum = 0, dn = 0;
    for (int b = 0; b < B; ++b) {
      double D = 0;
      for (int j = 0; j < T; ++j) {
        size_t o = (size_t)b * T + j;
        double x = ((double)lp[o] - (double)sp[o]) * (m[o] ? 1.0 : 0.0);
        x = std::isnan(x) ? 0.0 : (std::isinf(x) ? (x > 0 ? FLT_MAX_ : -FLT_MAX_) : x);
        x = (double)(float)x;
        if (m[o] && !std::isnan(sp[o])) { dsum += std::fabs((double)lp[o] - sp[o]); dn += 1; }
        if (seq) { D += x; continue; }
        double want = std::exp(x);
        want = (c.mode == TOKEN_TRUNCATE) ? std::min(std::max(want, (double)c.lo), (double)c.hi)
                                          : ((want < c.lo || want > c.hi) ? 0.0 : want);
        double got = r[o];
        if (!(std::fabs(got - want) <= 3e-7 * std::max(1.0, std::fabs(want)) || got == want)) ++bad;
      }
      if (seq) {
        double want = std::exp(D);
        want = (c.mode == SEQUENCE_TRUNCATE) ? std::min(std::max(want, (double)c.lo), (double)c.hi)
                                             : ((want < c.lo || want > c.hi) ? 0.0 : want);
        if (!(std::fabs(r[b] - want) <= 1e-4 * std::max(1.0, std::fabs(want)) || r[b] == want)) ++bad;
      }
    }
    CHECK(bad == 0, "is mode %d [%g, %g]: %d ratios differ", c.mode, c.lo, c.hi, bad);
    CHECK(st[1] == (float)dn, "d_n %g vs %g", st[1], dn);
    CHECK(st[0] == (float)dsum || std::fabs(st[0] - dsum) <= 1e-4 * dsum, "d_sum %g vs %g", st[0], dsum);
    if (!seq) {
      CHECK(r[(size_t)7 * T] == is_bound(1.0f, c.mode, c.lo, c.hi), "masked token ratio must be bound(1)");
    }
  }
  std::printf("is_ref: 6 mode/bound cases match the fp64 torch-semantics check\n");
}

int main() {
  test_exp_log();
  test_preprocess();
  test_rope();
  test_lse();
  test_is();
  std::printf(g_fail ? "test_cpu: %d FAILURES\n" : "test_cpu: all passed\n", g_fail);
  return g_fail ? 1 : 0;
}
