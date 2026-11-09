// K3 GPU API (launchers in k3_kernels.cu). All pointers are device pointers; every call is asynchronous on `stream`
// and makes no host synchronisation.
#pragma once

#include <cuda_runtime.h>
#include <stddef.h>
#include <stdint.h>

namespace k3gpu {

// ---------------- (a) image preprocessing ----------------
struct PreArgs {
  int H = 512, W = 512, C = 3;
  int ps = 16, ms = 2, tps = 2;
  bool hwc = true;                            // src (n_src, H, W, C) if true, else (n_src, C, H, W); uint8
  float mean[3] = {127.5f, 127.5f, 127.5f};   // HF's fused fp32 mean' (= mean / rescale_factor)
  float std_[3] = {127.5f, 127.5f, 127.5f};   // HF's fused fp32 std'
};
// out: (n_out * (H/ps) * (W/ps), C * tps * ps * ps), fp32 or bf16 (round to nearest even of the fp32 value).
// idx (n_out entries, may be null = identity): source image of each output image. Needs C == 3, ps % 4 == 0,
// W % 4 == 0, H and W multiples of ps * ms.
cudaError_t preprocess(const uint8_t* src, const int* idx, int n_out, const PreArgs& a, void* out, bool out_bf16,
                       cudaStream_t stream);

// ---------------- (c) rope positions ----------------
// pos: (3, B, W) int64. Template (computed once per prompt layout): ids, types (Lp each), pos (3 x Lp), end.
// err (one int, zeroed here): bit 1 = a row's prompt differs from the template, 2 = a completion token is not text,
// 4 = a row has fewer real tokens than the template. Positions are valid only where err == 0.
cudaError_t rope_positions(const int64_t* ids, const int64_t* mask, const int64_t* types, int B, int W,
                           const int64_t* tmpl_ids, const int64_t* tmpl_types, const int64_t* tmpl_pos, int64_t Lp,
                           int64_t tmpl_end, int64_t* pos, int* err, cudaStream_t stream);

// ---------------- (b) fused log-probabilities ----------------
// logp[n] = log softmax(inv_t * H[n] @ W^T)[tgt[n]] (0 where tgt is outside [0, V)), ent[n] the entropy (may be null).
// H (N, D) and W (V, D) bf16 (raw 16-bit storage), row-major, 16-byte aligned, D % 32 == 0. No (N, V) buffer:
// the workspace holds 4 * S * N floats. exact = true: CUDA-core fma chains (bit-equal to k3ref::lse_ref);
// exact = false: tensor cores (WMMA bf16, fp32 accumulate).
size_t lse_workspace_bytes(int N, int S);
cudaError_t lse_logprob(const uint16_t* H, const uint16_t* W, const int32_t* tgt, int N, int V, int D, float inv_t,
                        int S, bool exact, void* workspace, float* logp, float* ent, cudaStream_t stream);

// ---------------- (b) importance-sampling ratio ----------------
// lp, sp: (B, T) fp32 contiguous; mask: (B, T) uint8. mode: k3::IsMode. lo / hi: -inf / +inf for "None".
// ratio: (B, T) for token modes, (B) for sequence modes. stats: 7 floats (d_sum, d_n, d_max, r_sum, r_n, r_min,
// r_max). libdevice_exp: CUDA expf (bit-equal to torch.exp on CUDA) instead of k3_expf (bit-equal to the CPU ref).
size_t is_workspace_bytes(int B);
cudaError_t is_ratio(const float* lp, const float* sp, const uint8_t* mask, int B, int T, int mode, float lo, float hi,
                     bool libdevice_exp, float* ratio, void* workspace, float* stats, cudaStream_t stream);

}  // namespace k3gpu
