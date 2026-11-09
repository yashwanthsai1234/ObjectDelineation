// K3 CUDA kernels for the GRPO loop (Qwen3-VL-2B on one RTX 5090, sm_120). Build: nvcc -O3 -std=c++17
// -arch=sm_120 --fmad=false (see Makefile / GPU_RUN.md). --fmad=false keeps every a*b+c unfused unless written as
// fmaf, which is what makes the fp32 results bit-equal to the CPU reference; the tensor-core path is unaffected.
#include <cuda_bf16.h>
#include <mma.h>

#include <cub/block/block_scan.cuh>

#include "k3_common.cuh"
#include "k3_kernels.h"

using namespace k3;

namespace k3gpu {

// =====================================================================================================
// (a) preprocessing: one thread writes 4 consecutive features (same channel, temporal copy and patch row; 4 pixels)
// =====================================================================================================
template <bool HWC, bool BF16>
__global__ void __launch_bounds__(256) preprocess_kernel(const uint8_t* __restrict__ src, const int* __restrict__ idx,
                                                         int n_out, PatchGeom g, float m0, float m1, float m2,
                                                         float s0, float s1, float s2, void* __restrict__ out) {
  const int F4 = g.F / 4;
  const long long total = (long long)n_out * g.G * F4;
  const size_t img_bytes = (size_t)g.H * g.W * 3;
  for (long long i = blockIdx.x * (long long)blockDim.x + threadIdx.x; i < total;
       i += (long long)gridDim.x * blockDim.x) {
    const int f = (int)(i % F4) * 4;
    const long long rest = i / F4;
    const int p = (int)(rest % g.G);
    const int n = (int)(rest / g.G);
    int c, y, x;
    k3_patch_src(g, p, f, &c, &y, &x);
    const uint8_t* img = src + (size_t)(idx ? __ldg(idx + n) : n) * img_bytes;
    uint8_t v[4];
    if (HWC) {  // pixels x..x+3, all channels: 12 bytes, 4-byte aligned (x % 4 == 0)
      const uint32_t* w = reinterpret_cast<const uint32_t*>(img + ((size_t)y * g.W + x) * 3);
      const uint32_t a = __ldg(w), b = __ldg(w + 1), d = __ldg(w + 2);
      const unsigned long long lo = (unsigned long long)a | ((unsigned long long)b << 32);
#pragma unroll
      for (int k = 0; k < 4; ++k) {  // byte 3k + c of the 12 (registers only, no local array)
        const int bi = 3 * k + c;
        v[k] = bi < 8 ? (uint8_t)(lo >> (8 * bi)) : (uint8_t)(d >> (8 * (bi - 8)));
      }
    } else {
      uint32_t a = __ldg(reinterpret_cast<const uint32_t*>(img + ((size_t)c * g.H + y) * g.W + x));
#pragma unroll
      for (int k = 0; k < 4; ++k) v[k] = (a >> (8 * k)) & 0xff;
    }
    const float mean = c == 0 ? m0 : (c == 1 ? m1 : m2);
    const float sd = c == 0 ? s0 : (c == 1 ? s1 : s2);
    float r[4];
#pragma unroll
    for (int k = 0; k < 4; ++k) r[k] = k3_norm(v[k], mean, sd);
    const size_t o = ((size_t)n * g.G + p) * g.F + f;
    if (BF16) {
      uint2 pk;
      pk.x = (uint32_t)f32_to_bf16_rne(r[0]) | ((uint32_t)f32_to_bf16_rne(r[1]) << 16);
      pk.y = (uint32_t)f32_to_bf16_rne(r[2]) | ((uint32_t)f32_to_bf16_rne(r[3]) << 16);
      *reinterpret_cast<uint2*>(static_cast<uint16_t*>(out) + o) = pk;
    } else {
      *reinterpret_cast<float4*>(static_cast<float*>(out) + o) = make_float4(r[0], r[1], r[2], r[3]);
    }
  }
}

cudaError_t preprocess(const uint8_t* src, const int* idx, int n_out, const PreArgs& a, void* out, bool out_bf16,
                       cudaStream_t stream) {
  if (a.C != 3 || a.ps % 4 || a.W % 4 || a.H % (a.ps * a.ms) || a.W % (a.ps * a.ms)) return cudaErrorInvalidValue;
  PatchGeom g;
  g.H = a.H; g.W = a.W; g.ps = a.ps; g.ms = a.ms; g.tps = a.tps;
  g.gh = a.H / a.ps; g.gw = a.W / a.ps; g.F = a.C * a.tps * a.ps * a.ps; g.G = g.gh * g.gw;
  const long long total = (long long)n_out * g.G * (g.F / 4);
  if (total == 0) return cudaSuccess;
  int dev, sms;
  cudaGetDevice(&dev);
  cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, dev);
  const int threads = 256;
  const long long want = (total + threads - 1) / threads;
  const int blocks = (int)(want < (long long)sms * 16 ? want : (long long)sms * 16);
#define K3_PRE(HWC, BF) \
  preprocess_kernel<HWC, BF><<<blocks, threads, 0, stream>>>(src, idx, n_out, g, a.mean[0], a.mean[1], a.mean[2], \
                                                             a.std_[0], a.std_[1], a.std_[2], out)
  if (a.hwc) {
    if (out_bf16) K3_PRE(true, true); else K3_PRE(true, false);
  } else {
    if (out_bf16) K3_PRE(false, true); else K3_PRE(false, false);
  }
#undef K3_PRE
  return cudaGetLastError();
}

// =====================================================================================================
// (c) rope positions: one block per row, block-wide exclusive scan of the mask gives each real token's number r
// =====================================================================================================
constexpr int ROPE_THREADS = 256;

__global__ void __launch_bounds__(ROPE_THREADS) rope_kernel(const int64_t* __restrict__ ids,
                                                           const int64_t* __restrict__ mask,
                                                           const int64_t* __restrict__ types, int B, int W,
                                                           const int64_t* __restrict__ tids,
                                                           const int64_t* __restrict__ ttypes,
                                                           const int64_t* __restrict__ tpos, int64_t Lp,
                                                           int64_t tend, int64_t* __restrict__ pos, int* err) {
  using Scan = cub::BlockScan<int, ROPE_THREADS>;
  __shared__ typename Scan::TempStorage tmp;
  const int b = blockIdx.x;
  const size_t BW = (size_t)B * W;
  int64_t carry = 0;
  int e = 0;
  for (int j0 = 0; j0 < W; j0 += ROPE_THREADS) {
    const int j = j0 + threadIdx.x;
    const size_t o = (size_t)b * W + j;
    const int m = (j < W && mask[o] != 0) ? 1 : 0;
    int excl, total;
    Scan(tmp).ExclusiveSum(m, excl, total);
    if (j < W) {
      if (m) {
        const int64_t r = carry + excl;
        if (r < Lp) {
          if (ids[o] != tids[r] || types[o] != ttypes[r]) e |= 1;
        } else if (types[o] != 0) {
          e |= 2;
        }
#pragma unroll
        for (int a = 0; a < 3; ++a) pos[a * BW + o] = k3_rope_pos(r, Lp, tpos + a * Lp, tend);
      } else {
#pragma unroll
        for (int a = 0; a < 3; ++a) pos[a * BW + o] = 0;
      }
    }
    carry += total;
    __syncthreads();  // TempStorage reuse
  }
  if (threadIdx.x == 0 && carry < Lp) e |= 4;
  if (e) atomicOr(err, e);
}

cudaError_t rope_positions(const int64_t* ids, const int64_t* mask, const int64_t* types, int B, int W,
                           const int64_t* tmpl_ids, const int64_t* tmpl_types, const int64_t* tmpl_pos, int64_t Lp,
                           int64_t tmpl_end, int64_t* pos, int* err, cudaStream_t stream) {
  cudaError_t st = cudaMemsetAsync(err, 0, sizeof(int), stream);
  if (st != cudaSuccess || B == 0) return st;
  rope_kernel<<<B, ROPE_THREADS, 0, stream>>>(ids, mask, types, B, W, tmpl_ids, tmpl_types, tmpl_pos, Lp, tmpl_end,
                                              pos, err);
  return cudaGetLastError();
}

// =====================================================================================================
// (b) fused log-probabilities. CTA = 64 tokens x one vocab split; per 128-column vocab tile: GEMM tile into shared
// memory (tensor cores or exact fma chains), then the fixed-order online-logsumexp epilogue. Only per-split
// (max, sumexp, sum e*z, target logit) leave the SM.
// =====================================================================================================
namespace wmma = nvcuda::wmma;
constexpr int LDA = LSE_BK + 8;   // bf16 elements per staged row (80 bytes: 16-byte aligned, fewer bank conflicts)
constexpr int LDC = LSE_BN + 4;   // floats per staged accumulator row
constexpr int LSE_THREADS = 256;  // 8 warps: 2 (rows) x 4 (cols) of 32 x 32 warp tiles; 4 lanes per row in the epilogue
constexpr int AB_BYTES = 2 * LSE_BM * LDA * 2 + 2 * LSE_BN * LDA * 2;  // double-buffered A and B: 30,720
constexpr int C_BYTES = LSE_BM * LDC * 4;                               // 33,792
constexpr int LSE_SMEM = AB_BYTES > C_BYTES ? AB_BYTES : C_BYTES;       // C aliases A/B (static, < 48 KB)

__device__ __forceinline__ void cp_async16(void* smem, const void* gmem, bool pred) {
  const unsigned s = (unsigned)__cvta_generic_to_shared(smem);
  const int sz = pred ? 16 : 0;  // src-size 0: zero-fill, no global read
  asm volatile("cp.async.cg.shared.global [%0], [%1], 16, %2;\n" ::"r"(s), "l"(gmem), "r"(sz));
}
__device__ __forceinline__ void cp_async_commit() { asm volatile("cp.async.commit_group;\n" ::); }
template <int N>
__device__ __forceinline__ void cp_async_wait() { asm volatile("cp.async.wait_group %0;\n" ::"n"(N)); }

// Stage the (64 x 32) H block and the (128 x 32) W block of k-step k0 into buffer `buf`. Out-of-range rows: zeros.
__device__ __forceinline__ void lse_load(const __nv_bfloat16* H, const __nv_bfloat16* W, int N, int V, int D, int n0,
                                         int v0, int k0, __nv_bfloat16* As, __nv_bfloat16* Bs, int buf) {
  const int tid = threadIdx.x;
  {
    const int r = tid >> 2, ch = tid & 3;
    const bool ok = n0 + r < N;
    const __nv_bfloat16* g = ok ? H + (size_t)(n0 + r) * D + k0 + ch * 8 : H;
    cp_async16(As + (buf * LSE_BM + r) * LDA + ch * 8, g, ok);
  }
#pragma unroll
  for (int i = 0; i < 2; ++i) {
    const int id = tid + i * LSE_THREADS, r = id >> 2, ch = id & 3;
    const bool ok = v0 + r < V;
    const __nv_bfloat16* g = ok ? W + (size_t)(v0 + r) * D + k0 + ch * 8 : W;
    cp_async16(Bs + (buf * LSE_BN + r) * LDA + ch * 8, g, ok);
  }
}

// Fixed-order epilogue for one row (4 lanes): fold the tile's z = Cs * inv_t into the running state.
__device__ __forceinline__ void lse_epilogue(const float* Cs, int row, int q, int nv, float inv_t, LseState& st,
                                             int tgt_local, float& tgt_z) {
  const float* cr = Cs + row * LDC;
  float z[LSE_PER_LANE];
  float mq = -K3_INF;
#pragma unroll
  for (int i = 0; i < LSE_PER_LANE; ++i) {
    const int j = LSE_LANES * i + q;
    z[i] = j < nv ? K3_MUL(cr[j], inv_t) : -K3_INF;
    if (j < nv) mq = k3_max(mq, z[i]);
  }
  mq = k3_max(mq, __shfl_xor_sync(0xffffffffu, mq, 1));
  mq = k3_max(mq, __shfl_xor_sync(0xffffffffu, mq, 2));
  const float new_m = k3_max(st.m, mq);
  float E = 0.0f, X = 0.0f;
  if (new_m != -K3_INF) {
#pragma unroll
    for (int i = 0; i < LSE_PER_LANE; ++i)
      if (LSE_LANES * i + q < nv) lse_lane_accumulate(z[i], new_m, E, X);
  }
  E = K3_ADD(E, __shfl_xor_sync(0xffffffffu, E, 1));  // lane 0: e0 + e1 (lane 1: e1 + e0, the same bits)
  X = K3_ADD(X, __shfl_xor_sync(0xffffffffu, X, 1));
  E = K3_ADD(E, __shfl_xor_sync(0xffffffffu, E, 2));  // (e0 + e1) + (e2 + e3) on all four lanes
  X = K3_ADD(X, __shfl_xor_sync(0xffffffffu, X, 2));
  lse_fold(st, new_m, E, X);
  if (tgt_local >= 0 && tgt_local < nv) tgt_z = K3_MUL(cr[tgt_local], inv_t);
}

template <bool EXACT>
__global__ void __launch_bounds__(LSE_THREADS) lse_kernel(const __nv_bfloat16* __restrict__ H,
                                                          const __nv_bfloat16* __restrict__ W,
                                                          const int32_t* __restrict__ tgt, int N, int V, int D,
                                                          float inv_t, int TPS, float* __restrict__ pm,
                                                          float* __restrict__ psum, float* __restrict__ px,
                                                          float* __restrict__ pt) {
  __shared__ __align__(128) unsigned char smem[LSE_SMEM];
  __nv_bfloat16* As = reinterpret_cast<__nv_bfloat16*>(smem);
  __nv_bfloat16* Bs = As + 2 * LSE_BM * LDA;
  float* Cs = reinterpret_cast<float*>(smem);
  const int tid = threadIdx.x, warp = tid >> 5, wm = warp >> 2, wn = warp & 3;
  const int row = tid >> 2, q = tid & 3;
  const int n0 = blockIdx.x * LSE_BM, s = blockIdx.y, N_ = N;
  const int n = n0 + row;
  const int tn = n < N ? tgt[n] : -1;
  const int NT = (V + LSE_BN - 1) / LSE_BN;
  const int t_begin = s * TPS, t_end = min(NT, t_begin + TPS);
  const int KT = D / LSE_BK;
  LseState st = lse_init();
  float tgt_z = u2f(0x7fc00000u);

  for (int t = t_begin; t < t_end; ++t) {
    const int v0 = t * LSE_BN;
    wmma::fragment<wmma::accumulator, 16, 16, 16, float> acc[2][2];
    float ex[EXACT ? LSE_PER_LANE : 1];
    if (EXACT) {
#pragma unroll
      for (int i = 0; i < (EXACT ? LSE_PER_LANE : 1); ++i) ex[i] = 0.0f;
    } else {
#pragma unroll
      for (int i = 0; i < 2; ++i)
#pragma unroll
        for (int j = 0; j < 2; ++j) wmma::fill_fragment(acc[i][j], 0.0f);
    }
    lse_load(H, W, N, V, D, n0, v0, 0, As, Bs, 0);
    cp_async_commit();
    for (int kt = 0; kt < KT; ++kt) {
      const int buf = kt & 1;
      if (kt + 1 < KT) {
        lse_load(H, W, N, V, D, n0, v0, (kt + 1) * LSE_BK, As, Bs, buf ^ 1);
        cp_async_commit();
        cp_async_wait<1>();
      } else {
        cp_async_wait<0>();
      }
      __syncthreads();
      const __nv_bfloat16* a = As + buf * LSE_BM * LDA;
      const __nv_bfloat16* b = Bs + buf * LSE_BN * LDA;
      if (EXACT) {
        // thread (row, q): columns j = 4i + q, k in increasing order: acc = fma(h_k, w_k, acc), as k3ref::lse_ref
        const __nv_bfloat16* ar = a + row * LDA;
#pragma unroll 4
        for (int kk = 0; kk < LSE_BK; ++kk) {
          const float h = __bfloat162float(ar[kk]);
#pragma unroll
          for (int i = 0; i < (EXACT ? LSE_PER_LANE : 1); ++i)
            ex[i] = K3_FMA(h, __bfloat162float(b[(LSE_LANES * i + q) * LDA + kk]), ex[i]);
        }
      } else {
#pragma unroll
        for (int kk = 0; kk < LSE_BK; kk += 16) {
          wmma::fragment<wmma::matrix_a, 16, 16, 16, __nv_bfloat16, wmma::row_major> fa[2];
          wmma::fragment<wmma::matrix_b, 16, 16, 16, __nv_bfloat16, wmma::col_major> fb[2];
#pragma unroll
          for (int i = 0; i < 2; ++i) wmma::load_matrix_sync(fa[i], a + (wm * 32 + i * 16) * LDA + kk, LDA);
#pragma unroll
          for (int j = 0; j < 2; ++j) wmma::load_matrix_sync(fb[j], b + (wn * 32 + j * 16) * LDA + kk, LDA);
#pragma unroll
          for (int i = 0; i < 2; ++i)
#pragma unroll
            for (int j = 0; j < 2; ++j) wmma::mma_sync(acc[i][j], fa[i], fb[j], acc[i][j]);
        }
      }
      __syncthreads();
    }
    if (EXACT) {
#pragma unroll
      for (int i = 0; i < (EXACT ? LSE_PER_LANE : 1); ++i) Cs[row * LDC + LSE_LANES * i + q] = ex[i];
    } else {
#pragma unroll
      for (int i = 0; i < 2; ++i)
#pragma unroll
        for (int j = 0; j < 2; ++j)
          wmma::store_matrix_sync(Cs + (wm * 32 + i * 16) * LDC + wn * 32 + j * 16, acc[i][j], LDC,
                                  wmma::mem_row_major);
    }
    __syncthreads();
    lse_epilogue(Cs, row, q, min(LSE_BN, V - v0), inv_t, st, tn - v0, tgt_z);
    __syncthreads();  // the next tile's loads overwrite Cs
  }
  if (q == 0 && n < N_) {
    const size_t o = (size_t)s * N_ + n;
    pm[o] = st.m;
    psum[o] = st.s;
    px[o] = st.x;
    pt[o] = tgt_z;
  }
}

__global__ void lse_finish_kernel(const float* __restrict__ pm, const float* __restrict__ psum,
                                  const float* __restrict__ px, const float* __restrict__ pt,
                                  const int32_t* __restrict__ tgt, int N, int V, int S, int TPS,
                                  float* __restrict__ logp, float* __restrict__ ent) {
  const int n = blockIdx.x * blockDim.x + threadIdx.x;
  if (n >= N) return;
  const int tn = tgt[n];
  const bool valid = tn >= 0 && tn < V;
  const float tz = valid ? pt[(size_t)((tn / LSE_BN) / TPS) * N + n] : 0.0f;
  const LseOut o = lse_finish(pm + n, psum + n, px + n, S, N, tz, valid);
  logp[n] = o.logp;
  if (ent) ent[n] = o.entropy;
}

size_t lse_workspace_bytes(int N, int S) { return (size_t)4 * S * N * sizeof(float); }

cudaError_t lse_logprob(const uint16_t* H, const uint16_t* W, const int32_t* tgt, int N, int V, int D, float inv_t,
                        int S, bool exact, void* workspace, float* logp, float* ent, cudaStream_t stream) {
  if (D % LSE_BK || S < 1 || V < 1) return cudaErrorInvalidValue;
  if (((uintptr_t)H | (uintptr_t)W) & 15) return cudaErrorMisalignedAddress;
  if (N == 0) return cudaSuccess;
  const int NT = (V + LSE_BN - 1) / LSE_BN;
  const int TPS = (NT + S - 1) / S;
  float* pm = static_cast<float*>(workspace);
  float* psum = pm + (size_t)S * N;
  float* px = psum + (size_t)S * N;
  float* pt = px + (size_t)S * N;
  dim3 grid((N + LSE_BM - 1) / LSE_BM, S);
  const auto* h = reinterpret_cast<const __nv_bfloat16*>(H);
  const auto* w = reinterpret_cast<const __nv_bfloat16*>(W);
  if (exact)
    lse_kernel<true><<<grid, LSE_THREADS, 0, stream>>>(h, w, tgt, N, V, D, inv_t, TPS, pm, psum, px, pt);
  else
    lse_kernel<false><<<grid, LSE_THREADS, 0, stream>>>(h, w, tgt, N, V, D, inv_t, TPS, pm, psum, px, pt);
  cudaError_t e = cudaGetLastError();
  if (e != cudaSuccess) return e;
  lse_finish_kernel<<<(N + 255) / 256, 256, 0, stream>>>(pm, psum, px, pt, tgt, N, V, S, TPS, logp, ent);
  return cudaGetLastError();
}

// =====================================================================================================
// (b) importance-sampling ratio: one block per row; fixed-shape tree reductions (replayed by k3ref::is_ref)
// =====================================================================================================
__device__ __forceinline__ float block_tree(float v, int op, float* red) {
  red[threadIdx.x] = v;
  __syncthreads();
#pragma unroll
  for (int stride = IS_THREADS / 2; stride > 0; stride >>= 1) {
    if ((int)threadIdx.x < stride) red[threadIdx.x] = red_apply(op, red[threadIdx.x], red[threadIdx.x + stride]);
    __syncthreads();
  }
  const float r = red[0];
  __syncthreads();
  return r;
}

__constant__ int K3_STAT_OP[IS_NSTATS] = {RED_SUM, RED_SUM, RED_MAX, RED_SUM, RED_SUM, RED_MIN, RED_MAX};

template <bool LIBDEV>
__device__ __forceinline__ float is_exp(float x) {
  return LIBDEV ? expf(x) : k3_expf(x);
}

template <bool LIBDEV>
__global__ void __launch_bounds__(IS_THREADS) is_rows_kernel(const float* __restrict__ lp,
                                                             const float* __restrict__ sp,
                                                             const uint8_t* __restrict__ mask, int T, int mode,
                                                             float lo, float hi, float* __restrict__ ratio,
                                                             float* __restrict__ row_stats) {
  __shared__ float red[IS_THREADS];
  const int b = blockIdx.x;
  const bool seq = mode == SEQUENCE_TRUNCATE || mode == SEQUENCE_MASK;
  float a[IS_NSTATS] = {0.0f, 0.0f, -K3_INF, 0.0f, 0.0f, K3_INF, -K3_INF};
  float dacc = 0.0f;
  for (int j = threadIdx.x; j < T; j += IS_THREADS) {
    const size_t o = (size_t)b * T + j;
    const float l = lp[o], s = sp[o];
    const uint8_t m = mask[o];
    const float d = is_diff(l, s, m);
    if (m && !k3_isnan(s)) {
      float ad = K3_SUB(l, s);
      ad = ad < 0.0f ? -ad : ad;
      a[0] = K3_ADD(a[0], ad);
      a[1] = K3_ADD(a[1], 1.0f);
      a[2] = k3_max(a[2], ad);
    }
    if (seq) {
      dacc = K3_ADD(dacc, d);
    } else {
      const float r = is_bound(is_exp<LIBDEV>(d), mode, lo, hi);
      ratio[o] = r;
      if (m) {
        a[3] = K3_ADD(a[3], r);
        a[4] = K3_ADD(a[4], 1.0f);
        a[5] = k3_min(a[5], r);
        a[6] = k3_max(a[6], r);
      }
    }
  }
  float* rs = row_stats + (size_t)b * IS_NSTATS;
  float out[IS_NSTATS];
#pragma unroll
  for (int k = 0; k < 3; ++k) out[k] = block_tree(a[k], K3_STAT_OP[k], red);
  if (seq) {
    const float D = block_tree(dacc, RED_SUM, red);
    const float r = is_bound(is_exp<LIBDEV>(D), mode, lo, hi);
    if (threadIdx.x == 0) ratio[b] = r;
    out[3] = r; out[4] = 1.0f; out[5] = r; out[6] = r;
  } else {
#pragma unroll
    for (int k = 3; k < IS_NSTATS; ++k) out[k] = block_tree(a[k], K3_STAT_OP[k], red);
  }
  if (threadIdx.x == 0)
#pragma unroll
    for (int k = 0; k < IS_NSTATS; ++k) rs[k] = out[k];
}

__global__ void __launch_bounds__(IS_THREADS) is_stats_kernel(const float* __restrict__ row_stats, int B,
                                                              float* __restrict__ stats) {
  __shared__ float red[IS_THREADS];
  const float init[IS_NSTATS] = {0.0f, 0.0f, -K3_INF, 0.0f, 0.0f, K3_INF, -K3_INF};
#pragma unroll
  for (int k = 0; k < IS_NSTATS; ++k) {
    float a = init[k];
    for (int b = threadIdx.x; b < B; b += IS_THREADS) a = red_apply(K3_STAT_OP[k], a, row_stats[(size_t)b * IS_NSTATS + k]);
    const float r = block_tree(a, K3_STAT_OP[k], red);
    if (threadIdx.x == 0) stats[k] = r;
  }
}

size_t is_workspace_bytes(int B) { return (size_t)B * IS_NSTATS * sizeof(float); }

cudaError_t is_ratio(const float* lp, const float* sp, const uint8_t* mask, int B, int T, int mode, float lo, float hi,
                     bool libdevice_exp, float* ratio, void* workspace, float* stats, cudaStream_t stream) {
  if (mode < 0 || mode > 3) return cudaErrorInvalidValue;
  float* rows = static_cast<float*>(workspace);
  if (B > 0) {
    if (libdevice_exp)
      is_rows_kernel<true><<<B, IS_THREADS, 0, stream>>>(lp, sp, mask, T, mode, lo, hi, ratio, rows);
    else
      is_rows_kernel<false><<<B, IS_THREADS, 0, stream>>>(lp, sp, mask, T, mode, lo, hi, ratio, rows);
    cudaError_t e = cudaGetLastError();
    if (e != cudaSuccess) return e;
  }
  is_stats_kernel<<<1, IS_THREADS, 0, stream>>>(rows, B, stats);
  return cudaGetLastError();
}

}  // namespace k3gpu
