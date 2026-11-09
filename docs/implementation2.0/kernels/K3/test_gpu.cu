// K3 GPU correctness: every kernel against the CPU reference (byte equality where the design promises it, a
// measured error bound for the tensor-core log-probabilities). Build/run: make gpu && ./build/test_gpu
#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

#include "k3_kernels.h"
#include "k3_ref.h"

using namespace k3;
using namespace k3ref;

static int g_fail = 0;
#define CK(x)                                                                             \
  do {                                                                                    \
    cudaError_t e_ = (x);                                                                 \
    if (e_ != cudaSuccess) {                                                              \
      std::printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); \
      std::exit(2);                                                                       \
    }                                                                                     \
  } while (0)
#define CHECK(cond, ...)                                    \
  do {                                                      \
    if (!(cond)) {                                          \
      ++g_fail;                                             \
      std::printf("FAIL %s:%d: ", __FILE__, __LINE__);      \
      std::printf(__VA_ARGS__);                             \
      std::printf("\n");                                    \
    }                                                       \
  } while (0)

template <class T>
struct Dev {
  T* p = nullptr;
  size_t n = 0;
  explicit Dev(size_t n_) : n(n_) { CK(cudaMalloc(&p, std::max<size_t>(n, 1) * sizeof(T))); }
  Dev(const std::vector<T>& v) : Dev(v.size()) { CK(cudaMemcpy(p, v.data(), n * sizeof(T), cudaMemcpyHostToDevice)); }
  ~Dev() { cudaFree(p); }
  std::vector<T> get() const {
    std::vector<T> v(n);
    CK(cudaMemcpy(v.data(), p, n * sizeof(T), cudaMemcpyDeviceToHost));
    return v;
  }
};

static double secs(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

static void test_preprocess(std::mt19937& rng) {
  PreParams p;
  const PatchGeom g = geom(p);
  const size_t img = (size_t)p.H * p.W * 3;
  const int n_src = 36;
  std::vector<uint8_t> src(img * n_src);
  for (auto& v : src) v = (uint8_t)(rng() & 0xff);
  std::fill(src.begin(), src.begin() + img, 0);
  std::fill(src.begin() + img, src.begin() + 2 * img, 255);
  for (size_t i = 0; i < img; ++i) src[2 * img + i] = (uint8_t)((i * 7 + i / 3) & 0xff);
  for (size_t i = 0; i < img; ++i) src[3 * img + i] = (uint8_t)((((i / 3) % 512 / 16 + (i / 3) / 512 / 16) & 1) * 255);
  std::vector<int> idx(256);
  for (int r = 0; r < 256; ++r) idx[r] = (r / 8 + 5 * (r % 3)) % n_src;  // 8 rows per image, shuffled a little
  Dev<uint8_t> dsrc(src);
  Dev<int> didx(idx);
  for (int hwc = 1; hwc >= 0; --hwc) {
    p.hwc = hwc;
    k3gpu::PreArgs a;
    a.hwc = hwc;
    std::vector<uint8_t> s2 = src;
    if (!hwc) {  // the same images in CHW order
      for (int n = 0; n < n_src; ++n)
        for (int y = 0; y < p.H; ++y)
          for (int x = 0; x < p.W; ++x)
            for (int c = 0; c < 3; ++c)
              s2[n * img + ((size_t)c * p.H + y) * p.W + x] = src[n * img + ((size_t)y * p.W + x) * 3 + c];
    }
    Dev<uint8_t> ds(s2);
    const size_t out_n = (size_t)idx.size() * g.G * g.F;
    std::vector<float> ref(out_n);
    preprocess_ref(s2.data(), idx.data(), (int)idx.size(), p, ref.data());
    Dev<float> o32(out_n);
    Dev<uint16_t> o16(out_n);
    CK(k3gpu::preprocess(ds.p, didx.p, (int)idx.size(), a, o32.p, false, 0));
    CK(k3gpu::preprocess(ds.p, didx.p, (int)idx.size(), a, o16.p, true, 0));
    CK(cudaDeviceSynchronize());
    auto g32 = o32.get();
    auto g16 = o16.get();
    size_t bad32 = 0, bad16 = 0;
    for (size_t i = 0; i < out_n; ++i) {
      bad32 += std::memcmp(&g32[i], &ref[i], 4) != 0;
      bad16 += g16[i] != f32_to_bf16_rne(ref[i]);
    }
    CHECK(bad32 == 0 && bad16 == 0, "preprocess hwc=%d: %zu fp32 and %zu bf16 values differ", hwc, bad32, bad16);
    std::printf("(a) preprocess hwc=%d: 256 rows from %d images, fp32 %s, bf16 %s\n", hwc, n_src,
                bad32 ? "DIFFER" : "bit-equal", bad16 ? "DIFFER" : "bit-equal");
  }
}

static void test_rope(std::mt19937& rng) {
  const int merge = 2, pre = 4, post = 66, lh = 16, lw = 16;
  std::vector<int64_t> pids, ptypes;
  for (int i = 0; i < pre; ++i) { pids.push_back(1000 + i); ptypes.push_back(0); }
  for (int i = 0; i < lh * lw; ++i) { pids.push_back(151655); ptypes.push_back(1); }
  for (int i = 0; i < post; ++i) { pids.push_back(2000 + i); ptypes.push_back(0); }
  const int Lp = (int)pids.size();
  const Grid grid = {1, lh * merge, lw * merge};
  RopeTemplate t = rope_template(pids.data(), ptypes.data(), Lp, {grid}, merge);
  const int B = 256, P = Lp + 7, W = P + 1536;
  std::vector<int64_t> ids((size_t)B * W, 151643), mask(ids.size(), 0), types(ids.size(), 0);
  for (int b = 0; b < B; ++b) {
    int lead = (int)(rng() % 8), L = b == 0 ? 0 : (b == 1 ? 1536 : (int)(rng() % 1537));
    for (int j = 0; j < Lp; ++j) {
      size_t o = (size_t)b * W + lead + j;
      ids[o] = pids[j]; types[o] = ptypes[j]; mask[o] = 1;
    }
    for (int j = 0; j < L; ++j) {
      size_t o = (size_t)b * W + P + j;
      ids[o] = (int64_t)(rng() % 151643); mask[o] = 1;
    }
  }
  mask[(size_t)2 * W + P + 3] = 0;
  Dev<int64_t> di(ids), dm(mask), dt(types), dti(t.ids), dtt(t.types), dtp(t.pos);
  Dev<int64_t> dpos((size_t)3 * B * W);
  Dev<int> derr(1);
  auto run = [&](const Dev<int64_t>& I, const Dev<int64_t>& M, const Dev<int64_t>& T) {
    CK(k3gpu::rope_positions(I.p, M.p, T.p, B, W, dti.p, dtt.p, dtp.p, Lp, t.end, dpos.p, derr.p, 0));
    CK(cudaDeviceSynchronize());
    return derr.get()[0];
  };
  int err = run(di, dm, dt);
  auto gp = dpos.get();
  std::vector<int64_t> pe((size_t)3 * B * W), row((size_t)3 * W);
  int eerr = 0;
  rope_emul(B, W, ids.data(), mask.data(), types.data(), t, pe.data(), &eerr);
  size_t bad = 0, badref = 0;
  for (size_t i = 0; i < pe.size(); ++i) bad += gp[i] != pe[i];
  for (int b = 0; b < B; ++b) {
    rope_ref_row(ids.data() + (size_t)b * W, mask.data() + (size_t)b * W, types.data() + (size_t)b * W, W, {grid},
                 merge, row.data());
    for (int a = 0; a < 3; ++a)
      for (int j = 0; j < W; ++j) badref += gp[(size_t)a * B * W + (size_t)b * W + j] != row[(size_t)a * W + j];
  }
  CHECK(err == 0 && eerr == 0 && bad == 0 && badref == 0, "rope: err %d, %zu differ from emul, %zu from HF port", err,
        bad, badref);
  std::printf("(c) rope: %d rows x %d, Lp %d: GPU == host emulation (%zu differ) == HF-port reference (%zu differ)\n",
              B, W, Lp, bad, badref);
  auto i2 = ids; i2[(size_t)5 * W + P - 2] += 1;
  Dev<int64_t> di2(i2);
  CHECK(run(di2, dm, dt) & 1, "rope: prompt mismatch not flagged");
  auto t2 = types; t2[(size_t)1 * W + W - 1] = 1;
  Dev<int64_t> dt2(t2);
  CHECK(run(di, dm, dt2) & 2, "rope: image token in answer not flagged");
  auto m2 = mask; for (int j = 0; j < W; ++j) m2[(size_t)4 * W + j] = 0;
  Dev<int64_t> dm2(m2);
  CHECK(run(di, dm2, dt) & 4, "rope: short row not flagged");
}

static void fill_bf16(std::vector<uint16_t>& v, std::mt19937& rng, float scale) {
  std::normal_distribution<float> nd(0.0f, scale);
  for (auto& x : v) x = f32_to_bf16_rne(nd(rng));
}

static void lse_case(const char* name, int N, int V, int D, int S, float inv_t, std::mt19937& rng, bool tc_check) {
  std::vector<uint16_t> H((size_t)N * D), W((size_t)V * D);
  fill_bf16(H, rng, 1.0f);
  fill_bf16(W, rng, 0.05f);
  for (int k = 0; k < D; ++k) {
    W[(size_t)17 * D + k] = f32_to_bf16_rne(bf16_to_f32(H[(size_t)5 * D + k]) * 0.2f);  // row 5 peaked on 17
    H[(size_t)6 * D + k] = 0;                                                         // row 6 flat
  }
  std::vector<int32_t> tgt(N);
  for (int n = 0; n < N; ++n) tgt[n] = (int32_t)(rng() % V);
  tgt[0] = 0; tgt[1] = V - 1; tgt[2] = -100; tgt[3] = V; tgt[5] = 17;
  Dev<uint16_t> dH(H), dW(W);
  Dev<int32_t> dt(tgt);
  Dev<unsigned char> ws(k3gpu::lse_workspace_bytes(N, S));
  Dev<float> lp(N), en(N), lp2(N), en2(N);
  CK(k3gpu::lse_logprob(dH.p, dW.p, dt.p, N, V, D, inv_t, S, true, ws.p, lp.p, en.p, 0));
  CK(k3gpu::lse_logprob(dH.p, dW.p, dt.p, N, V, D, inv_t, S, false, ws.p, lp2.p, en2.p, 0));
  CK(cudaDeviceSynchronize());
  auto t0 = std::chrono::steady_clock::now();
  std::vector<float> rl(N), re(N);
  lse_ref(H.data(), W.data(), tgt.data(), N, V, D, inv_t, S, rl.data(), re.data());
  const double tref = secs(t0);
  auto gl = lp.get(), ge = en.get();
  size_t bad = 0;
  for (int n = 0; n < N; ++n) bad += std::memcmp(&gl[n], &rl[n], 4) != 0 || std::memcmp(&ge[n], &re[n], 4) != 0;
  CHECK(bad == 0, "lse %s exact: %zu of %d tokens differ from the CPU reference", name, bad, N);
  std::printf("(b) lse %-14s N %d V %d D %d S %d: exact mode %s the CPU reference (%zu differ; ref %.1f s)\n", name, N,
              V, D, S, bad ? "DIFFERS from" : "bit-equal to", bad, tref);
  if (tc_check) {
    std::vector<double> fl(N), fe(N);
    lse_ref_f64(H.data(), W.data(), tgt.data(), N, V, D, (double)inv_t, fl.data(), fe.data());
    auto tl = lp2.get(), te = en2.get();
    double dtl = 0, dte = 0, dxl = 0;
    for (int n = 0; n < N; ++n) {
      dtl = std::max(dtl, std::fabs(tl[n] - fl[n]));
      dte = std::max(dte, std::fabs(te[n] - fe[n]));
      dxl = std::max(dxl, std::fabs(gl[n] - fl[n]));
    }
    std::printf("    tensor-core mode vs fp64: max |dlogp| %.2e, max |dentropy| %.2e (exact mode vs fp64: %.2e)\n", dtl,
                dte, dxl);
    CHECK(dtl < 2e-3 && dte < 2e-3, "lse %s tensor-core error too large (%.2e, %.2e)", name, dtl, dte);
  }
}

static void test_is(std::mt19937& rng) {
  std::normal_distribution<float> nd(0.0f, 0.05f);
  const int B = 256, T = 1536;
  std::vector<float> lp((size_t)B * T), sp(lp.size());
  std::vector<uint8_t> m(lp.size());
  for (size_t i = 0; i < lp.size(); ++i) {
    lp[i] = -std::fabs(nd(rng)) * 20;
    sp[i] = lp[i] + nd(rng);
    m[i] = (rng() % 10) != 0;
  }
  sp[3] = NAN; m[4] = 0; sp[4] = NAN; sp[5] = -INFINITY; lp[6] = -80.0f; sp[6] = 0;
  for (int j = 0; j < T; ++j) m[(size_t)7 * T + j] = 0;
  Dev<float> dl(lp), ds(sp);
  Dev<uint8_t> dm(m);
  const float INF = K3_INF;
  struct Case { int mode; float lo, hi; } cases[] = {
      {TOKEN_TRUNCATE, -INF, 3.0f}, {TOKEN_MASK, 0.0f, 3.0f}, {TOKEN_TRUNCATE, 0.5f, 2.0f},
      {SEQUENCE_TRUNCATE, -INF, 3.0f}, {SEQUENCE_MASK, -INF, 3.0f}, {TOKEN_MASK, -INF, INF}};
  for (const Case& c : cases) {
    const bool seq = c.mode >= SEQUENCE_TRUNCATE;
    const size_t nr = seq ? B : (size_t)B * T;
    std::vector<float> rr(nr), rs((size_t)B * IS_NSTATS), st(IS_NSTATS);
    is_ref(lp.data(), sp.data(), m.data(), B, T, c.mode, c.lo, c.hi, rr.data(), rs.data(), st.data());
    Dev<float> dr(nr), dst(IS_NSTATS);
    Dev<unsigned char> ws(k3gpu::is_workspace_bytes(B));
    CK(k3gpu::is_ratio(dl.p, ds.p, dm.p, B, T, c.mode, c.lo, c.hi, false, dr.p, ws.p, dst.p, 0));
    CK(cudaDeviceSynchronize());
    auto gr = dr.get(), gs = dst.get();
    const bool ok = std::memcmp(gr.data(), rr.data(), nr * 4) == 0 && std::memcmp(gs.data(), st.data(), 28) == 0;
    CHECK(ok, "is mode %d [%g, %g]: GPU != CPU reference", c.mode, c.lo, c.hi);
    std::printf("(b) is_ratio mode %d [%g, %g]: ratios and 7 stats %s the CPU reference\n", c.mode, c.lo, c.hi,
                ok ? "bit-equal to" : "DIFFER from");
  }
}

int main(int argc, char** argv) {
  const bool quick = argc > 1 && std::strcmp(argv[1], "--quick") == 0;  // compute-sanitizer: skip the big case
  cudaDeviceProp prop;
  CK(cudaGetDeviceProperties(&prop, 0));
  std::printf("device: %s, cc %d.%d, %d SMs, %.1f GB\n", prop.name, prop.major, prop.minor, prop.multiProcessorCount,
              prop.totalGlobalMem / 1e9);
  std::mt19937 rng(7);
  test_preprocess(rng);
  test_rope(rng);
  lse_case("edge", 130, 1000, 64, 16, 1.0f, rng, true);   // 8 tiles, 16 splits (8 empty), N % 64 != 0
  lse_case("edge-T0.7", 77, 300, 96, 3, 1.0f / 0.7f, rng, true);
  if (!quick) lse_case("qwen3vl-2b", 192, 151936, 2048, 4, 1.0f, rng, true);  // the real head shape
  else lse_case("quick-2b", 64, 2048, 2048, 4, 1.0f, rng, true);
  test_is(rng);
  std::printf(g_fail ? "test_gpu: %d FAILURES\n" : "test_gpu: all passed\n", g_fail);
  return g_fail ? 1 : 0;
}
