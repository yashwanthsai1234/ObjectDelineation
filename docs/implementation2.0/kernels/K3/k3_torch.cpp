// PyTorch binding of the K3 kernels (torch.utils.cpp_extension.load; see test_torch_gpu.py). Every call runs on the
// current CUDA stream and reads nothing back to the host.
#include <ATen/cuda/CUDAContext.h>
#include <torch/extension.h>

#include <cmath>
#include <optional>

#include "k3_common.cuh"
#include "k3_kernels.h"

namespace {

void ok(cudaError_t e, const char* what) { TORCH_CHECK(e == cudaSuccess, what, ": ", cudaGetErrorString(e)); }
void dev(const torch::Tensor& t, const char* name) {
  TORCH_CHECK(t.is_cuda() && t.is_contiguous(), name, " must be a contiguous CUDA tensor");
}

torch::Tensor preprocess(torch::Tensor src, std::optional<torch::Tensor> idx, std::vector<double> mean,
                         std::vector<double> std_, bool hwc, bool bf16, int64_t ps, int64_t ms, int64_t tps) {
  dev(src, "src");
  TORCH_CHECK(src.scalar_type() == torch::kUInt8 && src.dim() == 4, "src: (n, H, W, 3) or (n, 3, H, W) uint8");
  TORCH_CHECK(mean.size() == 3 && std_.size() == 3, "mean/std: 3 values (HF's fused fp32 mean' and std')");
  k3gpu::PreArgs a;
  a.hwc = hwc;
  a.H = (int)(hwc ? src.size(1) : src.size(2));
  a.W = (int)(hwc ? src.size(2) : src.size(3));
  a.C = (int)(hwc ? src.size(3) : src.size(1));
  a.ps = (int)ps; a.ms = (int)ms; a.tps = (int)tps;
  for (int c = 0; c < 3; ++c) {
    a.mean[c] = (float)mean[c];
    a.std_[c] = (float)std_[c];
  }
  const int* ip = nullptr;
  int n_out = (int)src.size(0);
  if (idx) {
    dev(*idx, "idx");
    TORCH_CHECK(idx->scalar_type() == torch::kInt32, "idx: int32");
    ip = idx->data_ptr<int>();
    n_out = (int)idx->numel();
  }
  const int64_t G = (a.H / a.ps) * (a.W / a.ps), F = (int64_t)a.C * a.tps * a.ps * a.ps;
  auto out = torch::empty({n_out * G, F}, src.options().dtype(bf16 ? torch::kBFloat16 : torch::kFloat32));
  ok(k3gpu::preprocess(src.data_ptr<uint8_t>(), ip, n_out, a, out.data_ptr(), bf16,
                       at::cuda::getCurrentCUDAStream()), "k3 preprocess");
  return out;
}

std::vector<torch::Tensor> rope(torch::Tensor ids, torch::Tensor mask, torch::Tensor types, torch::Tensor tmpl_ids,
                                torch::Tensor tmpl_types, torch::Tensor tmpl_pos, int64_t tmpl_end) {
  for (auto* t : {&ids, &mask, &types, &tmpl_ids, &tmpl_types, &tmpl_pos}) {
    dev(*t, "rope inputs");
    TORCH_CHECK(t->scalar_type() == torch::kInt64, "rope inputs: int64");
  }
  const int B = (int)ids.size(0), W = (int)ids.size(1);
  const int64_t Lp = tmpl_ids.numel();
  TORCH_CHECK(tmpl_pos.numel() == 3 * Lp, "tmpl_pos: (3, Lp)");
  auto pos = torch::empty({3, B, W}, ids.options());
  auto err = torch::empty({1}, ids.options().dtype(torch::kInt32));
  ok(k3gpu::rope_positions(ids.data_ptr<int64_t>(), mask.data_ptr<int64_t>(), types.data_ptr<int64_t>(), B, W,
                           tmpl_ids.data_ptr<int64_t>(), tmpl_types.data_ptr<int64_t>(), tmpl_pos.data_ptr<int64_t>(),
                           Lp, tmpl_end, pos.data_ptr<int64_t>(), err.data_ptr<int>(),
                           at::cuda::getCurrentCUDAStream()), "k3 rope");
  return {pos, err};
}

int64_t default_splits(int64_t N) {
  int d, sms;
  cudaGetDevice(&d);
  cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, d);
  const int64_t tiles = (N + k3::LSE_BM - 1) / k3::LSE_BM;
  return std::min<int64_t>(64, std::max<int64_t>(1, (4 * sms + tiles - 1) / tiles));
}

std::vector<torch::Tensor> lse(torch::Tensor H, torch::Tensor W, torch::Tensor tgt, double inv_t, int64_t splits,
                               bool exact, bool entropy) {
  dev(H, "H"); dev(W, "W"); dev(tgt, "tgt");
  TORCH_CHECK(H.scalar_type() == torch::kBFloat16 && W.scalar_type() == torch::kBFloat16, "H, W: bf16");
  TORCH_CHECK(tgt.scalar_type() == torch::kInt32, "tgt: int32");
  const int N = (int)H.size(0), D = (int)H.size(1), V = (int)W.size(0);
  TORCH_CHECK(W.size(1) == D && tgt.numel() == N, "shapes");
  const int S = (int)(splits > 0 ? splits : default_splits(N));
  auto f = H.options().dtype(torch::kFloat32);
  auto ws = torch::empty({(int64_t)(k3gpu::lse_workspace_bytes(N, S) / 4)}, f);
  auto logp = torch::empty({N}, f);
  auto ent = entropy ? torch::empty({N}, f) : torch::Tensor();
  ok(k3gpu::lse_logprob(reinterpret_cast<const uint16_t*>(H.data_ptr()), reinterpret_cast<const uint16_t*>(W.data_ptr()),
                        tgt.data_ptr<int32_t>(), N, V, D, (float)inv_t, S, exact, ws.data_ptr(), logp.data_ptr<float>(),
                        entropy ? ent.data_ptr<float>() : nullptr, at::cuda::getCurrentCUDAStream()),
     "k3 lse");
  return {logp, ent};
}

std::vector<torch::Tensor> is_ratio(torch::Tensor lp, torch::Tensor sp, torch::Tensor mask, int64_t mode,
                                    std::optional<double> lo, std::optional<double> hi, bool libdevice) {
  dev(lp, "logps"); dev(sp, "sampling"); dev(mask, "mask");
  TORCH_CHECK(lp.scalar_type() == torch::kFloat32 && sp.scalar_type() == torch::kFloat32, "logps, sampling: fp32");
  TORCH_CHECK(mask.scalar_type() == torch::kBool || mask.scalar_type() == torch::kUInt8, "mask: bool or uint8");
  const int B = (int)lp.size(0), T = (int)lp.size(1);
  const bool seq = mode >= k3::SEQUENCE_TRUNCATE;
  auto f = lp.options();
  auto ratio = seq ? torch::empty({B, 1}, f) : torch::empty({B, T}, f);
  auto ws = torch::empty({(int64_t)(k3gpu::is_workspace_bytes(B) / 4)}, f);
  auto stats = torch::empty({k3::IS_NSTATS}, f);
  ok(k3gpu::is_ratio(lp.data_ptr<float>(), sp.data_ptr<float>(), reinterpret_cast<const uint8_t*>(mask.data_ptr()), B,
                     T, (int)mode, lo ? (float)*lo : -INFINITY, hi ? (float)*hi : INFINITY, libdevice,
                     ratio.data_ptr<float>(), ws.data_ptr(), stats.data_ptr<float>(),
                     at::cuda::getCurrentCUDAStream()),
     "k3 is_ratio");
  return {ratio, stats};
}

}  // namespace

PYBIND11_MODULE(TORCH_EXTENSION_NAME, m) {
  m.def("preprocess", &preprocess, "Qwen3-VL patches from uint8 images", py::arg("src"), py::arg("idx") = py::none(),
        py::arg("mean"), py::arg("std"), py::arg("hwc") = true, py::arg("bf16") = false, py::arg("ps") = 16,
        py::arg("ms") = 2, py::arg("tps") = 2);
  m.def("rope", &rope, "Qwen3-VL rope positions from a prompt template");
  m.def("lse", &lse, "fused per-token log-prob and entropy", py::arg("H"), py::arg("W"), py::arg("tgt"),
        py::arg("inv_t") = 1.0, py::arg("splits") = 0, py::arg("exact") = false, py::arg("entropy") = true);
  m.def("is_ratio", &is_ratio, "vLLM importance-sampling ratio (TRL modes) + its 7 stats", py::arg("logps"),
        py::arg("sampling"), py::arg("mask"), py::arg("mode"), py::arg("lo") = py::none(), py::arg("hi") = py::none(),
        py::arg("libdevice") = true);
  m.def("default_splits", &default_splits);
}
