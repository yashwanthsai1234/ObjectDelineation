// K3 micro-benchmarks (CUDA events, warm-up then median of reps). Shapes of one R5 GRPO step on Qwen3-VL-2B:
// 32 images x 8 answers = 256 rows, prompt ~325 tokens, answers up to 1,536 tokens, update micro-batch 16 rows.
// Build/run: make gpu && ./build/bench_gpu   (prints one line per case; ~1 minute on an RTX 5090)
#include <cuda_runtime.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

#include "k3_common.cuh"
#include "k3_kernels.h"

#define CK(x)                                                                             \
  do {                                                                                    \
    cudaError_t e_ = (x);                                                                 \
    if (e_ != cudaSuccess) {                                                              \
      std::printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); \
      std::exit(2);                                                                       \
    }                                                                                     \
  } while (0)

template <class F>
static float time_ms(F f, int reps) {
  cudaEvent_t a, b;
  CK(cudaEventCreate(&a));
  CK(cudaEventCreate(&b));
  for (int i = 0; i < 2; ++i) f();  // warm-up
  CK(cudaDeviceSynchronize());
  std::vector<float> t;
  for (int i = 0; i < reps; ++i) {
    CK(cudaEventRecord(a));
    f();
    CK(cudaEventRecord(b));
    CK(cudaEventSynchronize(b));
    float ms;
    CK(cudaEventElapsedTime(&ms, a, b));
    t.push_back(ms);
  }
  std::sort(t.begin(), t.end());
  return t[t.size() / 2];
}

static int default_splits(int N, int sms) {  // >= 4 CTAs per SM in flight, at most 64 splits
  const int tiles = (N + k3::LSE_BM - 1) / k3::LSE_BM;
  return std::min(64, std::max(1, (4 * sms + tiles - 1) / tiles));
}

int main() {
  cudaDeviceProp prop;
  CK(cudaGetDeviceProperties(&prop, 0));
  const int sms = prop.multiProcessorCount;
  std::printf("device: %s, cc %d.%d, %d SMs\n", prop.name, prop.major, prop.minor, sms);
  std::mt19937 rng(11);

  // ---------------- (a) preprocessing ----------------
  {
    const int n_src = 32, H = 512, W = 512;
    std::vector<uint8_t> src((size_t)n_src * H * W * 3);
    for (auto& v : src) v = (uint8_t)(rng() & 0xff);
    uint8_t* ds;
    CK(cudaMalloc(&ds, src.size()));
    CK(cudaMemcpy(ds, src.data(), src.size(), cudaMemcpyHostToDevice));
    std::vector<int> idx(256);
    for (int r = 0; r < 256; ++r) idx[r] = r / 8;
    int* di;
    CK(cudaMalloc(&di, idx.size() * 4));
    CK(cudaMemcpy(di, idx.data(), idx.size() * 4, cudaMemcpyHostToDevice));
    const size_t per = (size_t)1024 * 1536;
    void* out;
    CK(cudaMalloc(&out, 256 * per * 4));
    k3gpu::PreArgs a;
    struct C { const char* name; int n; const int* idx; bool bf16; } cs[] = {
        {"32 images -> fp32", 32, nullptr, false}, {"32 images -> bf16", 32, nullptr, true},
        {"256 rows  -> fp32", 256, di, false},     {"256 rows  -> bf16", 256, di, true}};
    for (const C& c : cs) {
      float ms = time_ms([&] { CK(k3gpu::preprocess(ds, c.idx, c.n, a, out, c.bf16, 0)); }, 20);
      const double bytes = (double)c.n * per * (c.bf16 ? 2 : 4) + (double)n_src * H * W * 3;
      std::printf("(a) preprocess %-18s %8.3f ms  %7.1f GB/s (output %.0f MB)\n", c.name, ms, bytes / ms / 1e6,
                  c.n * per * (c.bf16 ? 2 : 4) / 1e6);
    }
    // H2D of the 32 decoded images from pinned memory (what replaces the 1.61 GB fp32 host build)
    uint8_t* pin;
    CK(cudaMallocHost(&pin, src.size()));
    std::copy(src.begin(), src.end(), pin);
    float ms = time_ms([&] { CK(cudaMemcpyAsync(ds, pin, src.size(), cudaMemcpyHostToDevice, 0)); }, 20);
    std::printf("(a) H2D 32 uint8 images (pinned, %.1f MB)   %8.3f ms\n", src.size() / 1e6, ms);
    CK(cudaFreeHost(pin));
    CK(cudaFree(out));
    CK(cudaFree(di));
    CK(cudaFree(ds));
  }

  // ---------------- (c) rope positions ----------------
  {
    const int B = 256, Lp = 325, W = Lp + 7 + 1536;
    std::vector<int64_t> tids(Lp), ttypes(Lp, 0), tpos(3 * Lp);
    for (int i = 0; i < Lp; ++i) { tids[i] = 1000 + i; tpos[i] = tpos[Lp + i] = tpos[2 * Lp + i] = i; }
    std::vector<int64_t> ids((size_t)B * W, 0), mask(ids.size(), 0), types(ids.size(), 0);
    for (int b = 0; b < B; ++b) {
      int L = (int)(rng() % 1537);
      for (int j = 0; j < Lp; ++j) { ids[(size_t)b * W + 7 + j] = tids[j]; mask[(size_t)b * W + 7 + j] = 1; }
      for (int j = 0; j < L; ++j) mask[(size_t)b * W + Lp + 7 + j] = 1;
    }
    auto up = [](const std::vector<int64_t>& v) {
      int64_t* p;
      CK(cudaMalloc(&p, v.size() * 8));
      CK(cudaMemcpy(p, v.data(), v.size() * 8, cudaMemcpyHostToDevice));
      return p;
    };
    int64_t *di = up(ids), *dm = up(mask), *dt = up(types), *dti = up(tids), *dtt = up(ttypes), *dtp = up(tpos);
    int64_t* pos;
    int* err;
    CK(cudaMalloc(&pos, (size_t)3 * B * W * 8));
    CK(cudaMalloc(&err, 4));
    float ms = time_ms([&] { CK(k3gpu::rope_positions(di, dm, dt, B, W, dti, dtt, dtp, Lp, Lp, pos, err, 0)); }, 50);
    std::printf("(c) rope positions %d rows x %d          %8.3f ms\n", B, W, ms);
    for (int64_t* p : {di, dm, dt, dti, dtt, dtp, pos}) CK(cudaFree(p));
    CK(cudaFree(err));
  }

  // ---------------- (b) fused log-probabilities ----------------
  {
    const int V = 151936, D = 2048;
    const int Ns[] = {2048, 8192, 24576};  // 16 rows x {128, 512, 1536} answer tokens
    const int Nmax = 24576;
    std::vector<uint16_t> h((size_t)Nmax * D), w((size_t)V * D);
    std::normal_distribution<float> nd(0.0f, 1.0f);
    for (auto& x : h) x = k3::f32_to_bf16_rne(nd(rng));
    for (auto& x : w) x = k3::f32_to_bf16_rne(0.05f * nd(rng));
    std::vector<int32_t> tg(Nmax);
    for (auto& t : tg) t = (int32_t)(rng() % V);
    uint16_t *dh, *dw;
    int32_t* dt;
    float *lp, *en;
    CK(cudaMalloc(&dh, h.size() * 2));
    CK(cudaMalloc(&dw, w.size() * 2));
    CK(cudaMalloc(&dt, tg.size() * 4));
    CK(cudaMalloc(&lp, Nmax * 4));
    CK(cudaMalloc(&en, Nmax * 4));
    CK(cudaMemcpy(dh, h.data(), h.size() * 2, cudaMemcpyHostToDevice));
    CK(cudaMemcpy(dw, w.data(), w.size() * 2, cudaMemcpyHostToDevice));
    CK(cudaMemcpy(dt, tg.data(), tg.size() * 4, cudaMemcpyHostToDevice));
    for (int N : Ns) {
      const int S = default_splits(N, sms);
      void* ws;
      CK(cudaMalloc(&ws, k3gpu::lse_workspace_bytes(N, S)));
      float ms = time_ms([&] { CK(k3gpu::lse_logprob(dh, dw, dt, N, V, D, 1.0f, S, false, ws, lp, en, 0)); }, 5);
      const double tf = 2.0 * N * (double)V * D / (ms * 1e-3) / 1e12;
      std::printf("(b) lse tensor-core N %5d S %2d               %8.3f ms  %6.1f TFLOP/s  workspace %.2f MB\n", N, S,
                  ms, tf, k3gpu::lse_workspace_bytes(N, S) / 1e6);
      CK(cudaFree(ws));
    }
    {
      const int N = 256, S = default_splits(N, sms);
      void* ws;
      CK(cudaMalloc(&ws, k3gpu::lse_workspace_bytes(N, S)));
      float ms = time_ms([&] { CK(k3gpu::lse_logprob(dh, dw, dt, N, V, D, 1.0f, S, true, ws, lp, en, 0)); }, 3);
      std::printf("(b) lse exact (CUDA cores) N %5d S %2d        %8.3f ms  %6.2f TFLOP/s\n", N, S, ms,
                  2.0 * N * (double)V * D / (ms * 1e-3) / 1e12);
      CK(cudaFree(ws));
    }
    for (void* p : {(void*)dh, (void*)dw, (void*)dt, (void*)lp, (void*)en}) CK(cudaFree(p));
  }

  // ---------------- (b) IS ratio ----------------
  {
    const int B = 256, T = 1536;
    std::vector<float> l((size_t)B * T), s(l.size());
    std::vector<uint8_t> m(l.size(), 1);
    std::normal_distribution<float> nd(0.0f, 0.05f);
    for (size_t i = 0; i < l.size(); ++i) { l[i] = -std::abs(nd(rng)); s[i] = l[i] + nd(rng); }
    float *dl, *ds, *dr, *st;
    uint8_t* dm;
    void* ws;
    CK(cudaMalloc(&dl, l.size() * 4));
    CK(cudaMalloc(&ds, l.size() * 4));
    CK(cudaMalloc(&dr, l.size() * 4));
    CK(cudaMalloc(&dm, l.size()));
    CK(cudaMalloc(&st, 64));
    CK(cudaMalloc(&ws, k3gpu::is_workspace_bytes(B)));
    CK(cudaMemcpy(dl, l.data(), l.size() * 4, cudaMemcpyHostToDevice));
    CK(cudaMemcpy(ds, s.data(), s.size() * 4, cudaMemcpyHostToDevice));
    CK(cudaMemcpy(dm, m.data(), m.size(), cudaMemcpyHostToDevice));
    for (int mode : {0, 3}) {
      float ms = time_ms([&] { CK(k3gpu::is_ratio(dl, ds, dm, B, T, mode, -k3::K3_INF, 3.0f, true, dr, ws, st, 0)); },
                         50);
      std::printf("(b) is_ratio %s %d x %d              %8.3f ms\n", mode ? "sequence_mask " : "token_truncate", B, T, ms);
    }
    for (void* p : {(void*)dl, (void*)ds, (void*)dr, (void*)dm, (void*)st, ws}) CK(cudaFree(p));
  }
  return 0;
}
