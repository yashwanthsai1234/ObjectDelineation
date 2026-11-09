// K3 shared host/device code: one definition of every per-element formula, used by the CUDA kernels
// (k3_kernels.cu) and by the CPU reference and host emulation (k3_ref.cpp). Bit-equality between the GPU and the CPU
// rests on three rules:
//   1. every float operation is IEEE-754 binary32 round-to-nearest-even (no fast-math, no FTZ, no implicit FMA
//      contraction: nvcc --fmad=false and clang/gcc -ffp-contract=off; FMAs are written explicitly as fmaf);
//   2. exp and log are our own (k3_expf, k3_logf), built only from +, -, *, /, fmaf, rint and bit casts, so the host
//      and the device compute the same bits (CUDA's expf and libm's expf differ in the last bit);
//   3. every reduction runs in a fixed order that the CPU reference replays exactly.
// Integer work (patch layout, rope positions) is exact by construction.
#pragma once

#include <stdint.h>

#if defined(__CUDACC__)
#define K3_HD __host__ __device__ __forceinline__
#else
#define K3_HD inline
#endif

#if defined(__CUDA_ARCH__)
#define K3_ADD(a, b) __fadd_rn((a), (b))
#define K3_SUB(a, b) __fsub_rn((a), (b))
#define K3_MUL(a, b) __fmul_rn((a), (b))
#define K3_DIV(a, b) __fdiv_rn((a), (b))
#define K3_FMA(a, b, c) __fmaf_rn((a), (b), (c))
#define K3_RINT(a) rintf(a)
#else
#include <cmath>
#include <cstring>
#define K3_ADD(a, b) ((float)(a) + (float)(b))
#define K3_SUB(a, b) ((float)(a) - (float)(b))
#define K3_MUL(a, b) ((float)(a) * (float)(b))
#define K3_DIV(a, b) ((float)(a) / (float)(b))
#define K3_FMA(a, b, c) std::fmaf((a), (b), (c))
#define K3_RINT(a) std::nearbyintf(a)  // default rounding mode: half to even, like rintf on the device
#endif

namespace k3 {

constexpr float K3_INF = __builtin_huge_valf();
constexpr float FLT_MAX_ = 3.402823466e+38f;

K3_HD uint32_t f2u(float f) {
#if defined(__CUDA_ARCH__)
  return __float_as_uint(f);
#else
  uint32_t u;
  std::memcpy(&u, &f, 4);
  return u;
#endif
}
K3_HD float u2f(uint32_t u) {
#if defined(__CUDA_ARCH__)
  return __uint_as_float(u);
#else
  float f;
  std::memcpy(&f, &u, 4);
  return f;
#endif
}
K3_HD bool k3_isnan(float x) { return x != x; }
K3_HD float k3_max(float a, float b) { return a > b ? a : b; }  // same NaN behaviour on both sides
K3_HD float k3_min(float a, float b) { return a < b ? a : b; }

// bf16 <-> fp32. bf16 -> fp32 is exact (the high 16 bits). fp32 -> bf16 rounds to nearest even, as torch's .to().
K3_HD float bf16_to_f32(uint16_t b) { return u2f(((uint32_t)b) << 16); }
K3_HD uint16_t f32_to_bf16_rne(float f) {
  uint32_t u = f2u(f);
  if ((u & 0x7fffffffu) > 0x7f800000u) return (uint16_t)((u >> 16) | 0x40u);  // NaN stays NaN (quiet)
  uint32_t lsb = (u >> 16) & 1u;
  u += 0x7fffu + lsb;
  return (uint16_t)(u >> 16);
}

// 2^k for k in [-126, 127] as a normal float.
K3_HD float pow2i(int k) { return u2f((uint32_t)(k + 127) << 23); }

// exp(x) in binary32 from +, *, fmaf, rint only (max error measured by test_cpu: see DESIGN.md).
// Cody-Waite reduction x = n ln2 + r, |r| <= ln2/2, degree-7 Taylor polynomial in Horner form, scale by 2^n in two
// steps (so results near the overflow and subnormal limits are exact scalings).
K3_HD float k3_expf(float x) {
  if (k3_isnan(x)) return x;
  if (x > 88.72283935546875f) return K3_INF;
  if (x < -103.972084045410156f) return 0.0f;
  const float log2e = 1.44269502162933349609f;
  const float ln2_hi = 0.693145751953125f;        // 16 significant bits: n * ln2_hi exact for |n| <= 2^8
  const float ln2_lo = 1.42860676533018704e-06f;  // ln2 - ln2_hi
  float n = K3_RINT(K3_MUL(x, log2e));
  float r = K3_FMA(-n, ln2_hi, x);
  r = K3_FMA(-n, ln2_lo, r);
  float p = 1.98412698412698413e-04f;              // 1/5040
  p = K3_FMA(p, r, 1.38888888888888889e-03f);      // 1/720
  p = K3_FMA(p, r, 8.33333333333333333e-03f);      // 1/120
  p = K3_FMA(p, r, 4.16666666666666667e-02f);      // 1/24
  p = K3_FMA(p, r, 1.66666666666666667e-01f);      // 1/6
  p = K3_FMA(p, r, 0.5f);
  p = K3_FMA(p, r, 1.0f);
  p = K3_FMA(p, r, 1.0f);
  int ni = (int)n;
  int n1 = ni / 2;
  int n2 = ni - n1;
  return K3_MUL(K3_MUL(p, pow2i(n1)), pow2i(n2));
}

// log(x) in binary32: x = 2^e m, m in [sqrt(1/2), sqrt(2)), log m = 2 atanh(s), s = (m-1)/(m+1), odd series to s^11.
K3_HD float k3_logf(float x) {
  if (k3_isnan(x) || x < 0.0f) return u2f(0x7fc00000u);
  if (x == 0.0f) return -K3_INF;
  if (x == K3_INF) return K3_INF;
  int e = 0;
  if (x < 1.17549435e-38f) {  // subnormal: scale into the normal range (exact)
    x = K3_MUL(x, 8388608.0f);
    e = -23;
  }
  uint32_t u = f2u(x);
  e += (int)((u >> 23) & 0xffu) - 127;
  float m = u2f((u & 0x007fffffu) | 0x3f800000u);  // [1, 2)
  if (m > 1.41421356237309505f) {
    m = K3_MUL(m, 0.5f);
    e += 1;
  }
  float f = K3_SUB(m, 1.0f);  // exact (Sterbenz)
  float s = K3_DIV(f, K3_ADD(2.0f, f));
  float s2 = K3_MUL(s, s);
  float q = 0.181818181818181818f;  // 2/11
  q = K3_FMA(q, s2, 0.222222222222222222f);  // 2/9
  q = K3_FMA(q, s2, 0.285714285714285714f);  // 2/7
  q = K3_FMA(q, s2, 0.4f);                   // 2/5
  q = K3_FMA(q, s2, 0.666666666666666667f);  // 2/3
  // log m = 2s + s*s2*q, with 2s = f - s*f (more accurate than doubling the rounded s)
  float two_s = K3_FMA(-s, f, f);
  float lm = K3_FMA(K3_MUL(s, s2), q, two_s);
  const float ln2_hi = 0.693145751953125f;
  const float ln2_lo = 1.42860676533018704e-06f;
  float fe = (float)e;
  return K3_FMA(fe, ln2_hi, K3_FMA(fe, ln2_lo, lm));
}

// ---------------------------------------------------------------------------------------------------------------
// (a) Qwen3-VL image preprocessing. HF Qwen2VLImageProcessor (torchvision backend), for an image whose smart_resize
// size equals its own size (512 x 512, factor 32), computes in fp32:
//     y = (float(u8) - mean') / std',  mean' = fp32(mean) * fp32(1/rescale), std' likewise
// (TorchvisionBackend._fuse_mean_std_and_rescale_factor, then torchvision normalize: sub_ then div_ by tensors).
// mean' and std' are passed in as the fp32 values HF computes (127.5 for Qwen3-VL), so the device repeats the
// same two correctly rounded operations.
K3_HD float k3_norm(uint8_t v, float mean_f, float std_f) { return K3_DIV(K3_SUB((float)v, mean_f), std_f); }

struct PatchGeom {
  int H, W;          // image size (already the smart_resize size)
  int ps, ms, tps;   // patch_size 16, merge_size 2, temporal_patch_size 2
  int gh, gw;        // H / ps, W / ps
  int F;             // C * tps * ps * ps (1536)
  int G;             // gh * gw (1024)
};

// Output element (patch p, feature f) of one image -> (channel, y, x) of the source pixel. Layout of HF patchify:
// p = ((bh * (gw/ms) + bw) * ms + mh) * ms + mw,  f = ((c * tps + t) * ps + py) * ps + px,
// y = (bh * ms + mh) * ps + py, x = (bw * ms + mw) * ps + px; t is the temporal copy (same pixel).
K3_HD void k3_patch_src(const PatchGeom& g, int p, int f, int* c, int* y, int* x) {
  int px = f % g.ps;
  int py = (f / g.ps) % g.ps;
  int ch = f / (g.ps * g.ps * g.tps);
  int mw = p % g.ms;
  int mh = (p / g.ms) % g.ms;
  int bw = (p / (g.ms * g.ms)) % (g.gw / g.ms);
  int bh = p / (g.ms * g.ms * (g.gw / g.ms));
  *c = ch;
  *y = (bh * g.ms + mh) * g.ps + py;
  *x = (bw * g.ms + mw) * g.ps + px;
}

// ---------------------------------------------------------------------------------------------------------------
// (c) rope positions. For a row whose real tokens (attention_mask != 0, in column order) are the template prompt
// (Lp tokens) followed by text tokens, Qwen3-VL get_rope_index gives real token number r the position
//     tmpl_pos[a][r]            for r < Lp,
//     tmpl_end + (r - Lp)       for r >= Lp (text continues the last group; all three axes equal),
// and 0 at masked columns.
K3_HD int64_t k3_rope_pos(int64_t r, int64_t Lp, const int64_t* tmpl_pos_a, int64_t tmpl_end) {
  return r < Lp ? tmpl_pos_a[r] : tmpl_end + (r - Lp);
}

// ---------------------------------------------------------------------------------------------------------------
// (b) fused log-probabilities. Fixed reduction structure (the CPU reference replays it):
//   vocab tiles of LSE_BN columns; split s owns tiles [s*TPS, min(NT, (s+1)*TPS)); in a tile, lane q in 0..3 of a
//   row owns columns j = 4*i + q (i = 0..31), summed in i order; lanes combine as (l0 + l1) + (l2 + l3).
constexpr int LSE_BM = 64;
constexpr int LSE_BN = 128;
constexpr int LSE_BK = 32;
constexpr int LSE_LANES = 4;
constexpr int LSE_PER_LANE = LSE_BN / LSE_LANES;  // 32

struct LseState {
  float m, s, x;  // running max, sum exp(z - m), sum exp(z - m) * z
};

K3_HD LseState lse_init() { return LseState{-K3_INF, 0.0f, 0.0f}; }

// Fold one tile's lane sums into the running state. new_m = max(state.m, tile max) is computed by the caller (the
// lane partial sums E, X are relative to new_m).
K3_HD void lse_fold(LseState& st, float new_m, float E, float X) {
  if (new_m == -K3_INF) return;  // no valid column yet
  float sc = (st.m == -K3_INF) ? 0.0f : k3_expf(K3_SUB(st.m, new_m));
  st.s = K3_FMA(st.s, sc, E);
  st.x = K3_FMA(st.x, sc, X);
  st.m = new_m;
}

// One lane's partial sums for a tile: z values z_i (i = 0..31, valid flags), relative to new_m.
K3_HD void lse_lane_accumulate(float z, float new_m, float& E, float& X) {
  float e = k3_expf(K3_SUB(z, new_m));
  E = K3_ADD(E, e);
  X = K3_FMA(e, z, X);
}

// Combine the splits' states (fixed order s = 0..S-1) and the target logit into (logp, entropy).
struct LseOut {
  float logp, entropy;
};
K3_HD LseOut lse_finish(const float* m, const float* s, const float* x, int S, int stride, float tgt, bool valid) {
  float M = -K3_INF;
  for (int k = 0; k < S; ++k) M = k3_max(M, m[k * stride]);
  float Ssum = 0.0f, Xsum = 0.0f;
  for (int k = 0; k < S; ++k) {
    float mk = m[k * stride];
    float sc = (mk == -K3_INF) ? 0.0f : k3_expf(K3_SUB(mk, M));
    Ssum = K3_FMA(s[k * stride], sc, Ssum);
    Xsum = K3_FMA(x[k * stride], sc, Xsum);
  }
  float lg = k3_logf(Ssum);
  LseOut o;
  o.logp = valid ? K3_SUB(K3_SUB(tgt, M), lg) : 0.0f;
  o.entropy = K3_SUB(K3_ADD(M, lg), K3_DIV(Xsum, Ssum));
  return o;
}

// ---------------------------------------------------------------------------------------------------------------
// (b) importance-sampling ratio, mirroring objdel.train.fast.is_ratio (TRL 1.17 grpo_trainer 2776-2824):
//   diff = nan_to_num((logp - sampling) * mask, nan=0)  (+inf -> FLT_MAX, -inf -> -FLT_MAX, as torch.nan_to_num)
//   token_*: r = exp(diff); sequence_*: r = exp(sum_t diff) per row
//   *_truncate: clamp(r, lo, hi) (token_truncate = token-level TIS);  *_mask: r if lo <= r <= hi else 0
enum IsMode : int { TOKEN_TRUNCATE = 0, TOKEN_MASK = 1, SEQUENCE_TRUNCATE = 2, SEQUENCE_MASK = 3 };
constexpr int IS_THREADS = 256;
constexpr int IS_NSTATS = 7;  // d_sum, d_n, d_max, r_sum, r_n, r_min, r_max (fast.FastGRPO.set_is_ratio order)

K3_HD float nan_to_num(float v) {
  if (k3_isnan(v)) return 0.0f;
  if (v == K3_INF) return FLT_MAX_;
  if (v == -K3_INF) return -FLT_MAX_;
  return v;
}
K3_HD float is_diff(float lp, float sp, uint8_t mask) {
  return nan_to_num(K3_MUL(K3_SUB(lp, sp), mask ? 1.0f : 0.0f));
}
K3_HD float is_bound(float r, int mode, float lo, float hi) {
  if (mode == TOKEN_TRUNCATE || mode == SEQUENCE_TRUNCATE) {
    // torch.clamp(r, min=lo, max=hi) = min(max(r, lo), hi); NaN propagates
    if (k3_isnan(r)) return r;
    if (r < lo) r = lo;
    if (r > hi) r = hi;
    return r;
  }
  return (r < lo || r > hi) ? 0.0f : r;
}

// Fixed-shape tree reduction over IS_THREADS values (stride 128, 64, ..., 1), replayed by the CPU reference.
enum RedOp : int { RED_SUM = 0, RED_MAX = 1, RED_MIN = 2 };
K3_HD float red_apply(int op, float a, float b) {
  return op == RED_SUM ? K3_ADD(a, b) : (op == RED_MAX ? k3_max(a, b) : k3_min(a, b));
}

}  // namespace k3
