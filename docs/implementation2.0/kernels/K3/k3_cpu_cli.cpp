// K3 CPU command line for the Python harness (test_hf_cpu.py): runs the CPU reference and the host emulation of the
// kernels on raw binary files, so their outputs can be compared byte for byte with Hugging Face's.
//   k3_cpu_cli preprocess IN.u8 N H W HWC M0 M1 M2 S0 S1 S2 OUT.f32 OUT.bf16
//   k3_cpu_cli rope IN.i64 OUT.i64
// rope IN: B W Lp merge gt gh gw, tmpl_ids[Lp], tmpl_types[Lp], ids[B*W], mask[B*W], types[B*W] (all int64)
// rope OUT: err, pos_emul[3*B*W], pos_ref[3*B*W] (int64)
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include "k3_ref.h"

using namespace k3ref;

template <class T>
static std::vector<T> read_all(const char* path) {
  FILE* f = std::fopen(path, "rb");
  if (!f) { std::perror(path); std::exit(2); }
  std::fseek(f, 0, SEEK_END);
  long n = std::ftell(f);
  std::fseek(f, 0, SEEK_SET);
  std::vector<T> v(n / sizeof(T));
  if (std::fread(v.data(), sizeof(T), v.size(), f) != v.size()) { std::perror(path); std::exit(2); }
  std::fclose(f);
  return v;
}
template <class T>
static void write_all(const char* path, const T* p, size_t n) {
  FILE* f = std::fopen(path, "wb");
  if (!f || std::fwrite(p, sizeof(T), n, f) != n) { std::perror(path); std::exit(2); }
  std::fclose(f);
}

int main(int argc, char** argv) {
  if (argc >= 15 && std::string(argv[1]) == "preprocess") {
    auto src = read_all<uint8_t>(argv[2]);
    PreParams p;
    const int n = std::atoi(argv[3]);
    p.H = std::atoi(argv[4]);
    p.W = std::atoi(argv[5]);
    p.hwc = std::atoi(argv[6]) != 0;
    for (int c = 0; c < 3; ++c) {
      p.mean[c] = std::strtof(argv[7 + c], nullptr);
      p.std_[c] = std::strtof(argv[10 + c], nullptr);
    }
    const k3::PatchGeom g = geom(p);
    if (src.size() != (size_t)n * p.H * p.W * 3) { std::fprintf(stderr, "bad input size\n"); return 2; }
    std::vector<float> a((size_t)n * g.G * g.F), b(a.size());
    preprocess_ref(src.data(), nullptr, n, p, a.data());
    preprocess_emul(src.data(), nullptr, n, p, b.data());
    if (std::memcmp(a.data(), b.data(), a.size() * 4) != 0) { std::fprintf(stderr, "emul != ref\n"); return 1; }
    std::vector<uint16_t> h(b.size());
    for (size_t i = 0; i < b.size(); ++i) h[i] = k3::f32_to_bf16_rne(b[i]);
    write_all(argv[13], b.data(), b.size());
    write_all(argv[14], h.data(), h.size());
    return 0;
  }
  if (argc == 4 && std::string(argv[1]) == "rope") {
    auto in = read_all<int64_t>(argv[2]);
    const int64_t B = in[0], W = in[1], Lp = in[2], merge = in[3];
    const Grid grid = {in[4], in[5], in[6]};
    const int64_t* tids = in.data() + 7;
    const int64_t* ttypes = tids + Lp;
    const int64_t* ids = ttypes + Lp;
    const int64_t* mask = ids + B * W;
    const int64_t* types = mask + B * W;
    if ((int64_t)in.size() != 7 + 2 * Lp + 3 * B * W) { std::fprintf(stderr, "bad rope input size\n"); return 2; }
    RopeTemplate t = rope_template(tids, ttypes, (int)Lp, {grid}, (int)merge);
    std::vector<int64_t> out(1 + 6 * B * W);
    int err = 0;
    rope_emul((int)B, (int)W, ids, mask, types, t, out.data() + 1, &err);
    out[0] = err;
    std::vector<int64_t> row(3 * W);
    int64_t* ref = out.data() + 1 + 3 * B * W;
    for (int64_t b = 0; b < B; ++b) {
      try {
        rope_ref_row(ids + b * W, mask + b * W, types + b * W, (int)W, {grid}, (int)merge, row.data());
      } catch (const std::exception&) {  // a row HF itself cannot position (e.g. an image token in an answer)
        std::fill(row.begin(), row.end(), -1);
      }
      for (int a = 0; a < 3; ++a)
        std::memcpy(ref + a * B * W + b * W, row.data() + a * W, W * sizeof(int64_t));
    }
    write_all(argv[3], out.data(), out.size());
    return 0;
  }
  std::fprintf(stderr, "usage: see the header of k3_cpu_cli.cpp\n");
  return 2;
}
