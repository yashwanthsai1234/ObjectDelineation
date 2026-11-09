// K3 CPU reference. Plain C++17, compiles with clang++ on macOS and g++ on Linux (no CUDA needed).
#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "k3_common.cuh"

namespace k3ref {

// ---------------- (a) preprocessing ----------------
struct PreParams {
  int H = 512, W = 512, C = 3;
  int ps = 16, ms = 2, tps = 2;
  bool hwc = true;  // source layout: (n, H, W, C) uint8 (PIL / numpy decode) or (n, C, H, W)
  float mean[3] = {127.5f, 127.5f, 127.5f};  // HF's fused fp32 mean' = mean / rescale_factor
  float std_[3] = {127.5f, 127.5f, 127.5f};  // HF's fused fp32 std'
};
k3::PatchGeom geom(const PreParams& p);

// HF-structured reference: normalise the whole CHW image, then walk HF patchify's reshape/permute/expand as an
// explicit loop nest in output order. Independent of k3_patch_src. out: (n_out * G, F) fp32. idx (n_out entries,
// nullable = identity) picks the source image of each output image (row duplication of a generation batch).
void preprocess_ref(const uint8_t* src, const int* idx, int n_out, const PreParams& p, float* out);
// Host emulation of the GPU kernel's element mapping (k3_patch_src + k3_norm), for the CPU tests.
void preprocess_emul(const uint8_t* src, const int* idx, int n_out, const PreParams& p, float* out);

// ---------------- (c) rope positions ----------------
using Grid = std::array<int64_t, 3>;
// Port of transformers 5.20.1 Qwen3VLModel.get_rope_index for one row (image grids only). pos: 3 x W (row-major,
// axis stride W); masked columns get 0. Returns the row's mrope delta (max position + 1 - real token count).
int64_t rope_ref_row(const int64_t* ids, const int64_t* mask, const int64_t* types, int W,
                     const std::vector<Grid>& grids, int merge, int64_t* pos);

struct RopeTemplate {
  int64_t Lp = 0;
  std::vector<int64_t> ids, types;  // Lp each
  std::vector<int64_t> pos;         // 3 x Lp
  int64_t end = 0;                  // position after the prompt (current_pos)
};
// Template from one row's real prompt tokens (no padding): computed once per prompt layout.
RopeTemplate rope_template(const int64_t* ids, const int64_t* types, int Lp, const std::vector<Grid>& grids,
                           int merge);
// Host emulation of the GPU kernel: pos (3, B, W); err bit 1 = prompt mismatch, 2 = non-text completion token,
// 4 = fewer real tokens than the template.
void rope_emul(int B, int W, const int64_t* ids, const int64_t* mask, const int64_t* types, const RopeTemplate& t,
               int64_t* pos, int* err);

// ---------------- (b) fused log-probabilities ----------------
// Exact-mode reference: z = fma chain over k = 0..D-1 of bf16 products (fp32), times inv_t, then the fixed
// tile / lane / split reduction of k3_common.cuh. H: (N, D) bf16 bits, Wv: (V, D) bf16 bits. Targets outside
// [0, V) give logp 0. Multithreaded over rows (results do not depend on the thread count).
void lse_ref(const uint16_t* H, const uint16_t* Wv, const int32_t* tgt, int N, int V, int D, float inv_t, int S,
             float* logp, float* ent);
// fp64 reference (accuracy yardstick for the tensor-core mode).
void lse_ref_f64(const uint16_t* H, const uint16_t* Wv, const int32_t* tgt, int N, int V, int D, double inv_t,
                 double* logp, double* ent);
int lse_tiles(int V);
int lse_tiles_per_split(int V, int S);

// ---------------- (b) importance-sampling ratio ----------------
// Replays the GPU kernel: per row, IS_THREADS strided partials + fixed tree. ratio: (B, T) for token modes, (B) for
// sequence modes. row_stats: (B, 7). stats: (7), reduced over rows the same way.
void is_ref(const float* lp, const float* sp, const uint8_t* mask, int B, int T, int mode, float lo, float hi,
            float* ratio, float* row_stats, float* stats);

}  // namespace k3ref
