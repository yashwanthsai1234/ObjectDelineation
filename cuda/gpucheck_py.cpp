// nanobind binding of cuda/gpucheck.h as objdel._gpucheck (CPU build by default; `make cuda` builds it with nvcc).
// Inputs are the packed arrays of data/checker_eval (experiments/cuda_checker/prepare_data.py); they are checked here
// with polycore's caps before the GIL is released, so bad input raises ValueError and never reaches the kernels. The
// GIL is released during the checker call (other Python threads run; gpucheck's mutex orders concurrent calls) and
// taken back before the result arrays are made. Pattern and references: src/objdel/_polycore.cpp.
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/optional.h>

#include <optional>
#include <stdexcept>
#include <vector>

#include "gpucheck.h"
#include "polycore/polycore.hpp"

namespace nb = nanobind;
using gpucheck::i32;
using gpucheck::i64;
using gpucheck::u8;

template <class T, int D>
using In = nb::ndarray<const T, nb::ndim<D>, nb::c_contig, nb::device::cpu>;

namespace {

template <class T>
nb::ndarray<nb::numpy, T> to_numpy(std::vector<T>&& v, std::initializer_list<size_t> shape) {
  auto* p = new std::vector<T>(std::move(v));
  nb::capsule owner(p, [](void* q) noexcept { delete static_cast<std::vector<T>*>(q); });
  return nb::ndarray<nb::numpy, T>(p->data(), shape, owner);
}

// start[0] == 0, non-decreasing, start[k] == total, every segment at most `cap` long.
void check_offsets(const In<i64, 1>& start, i64 total, i64 cap, const char* what) {
  const i64 k = i64(start.shape(0)) - 1;
  if (k < 0 || start(0) != 0 || start(k) != total)
    throw std::invalid_argument(std::string(what) + " must run from 0 to the number of items it indexes");
  for (i64 i = 0; i < k; ++i)
    if (start(i + 1) < start(i) || start(i + 1) - start(i) > cap)
      throw std::invalid_argument(std::string(what) + " must be non-decreasing, within the caps");
}

}  // namespace

NB_MODULE(_gpucheck, m) {
  m.doc() = "Batched exact checker (cuda/gpucheck.h): polycore's even-odd masks (bit-packed rows) and validate "
            "totals, polygonize's rings.";
  m.def("build", [] { return gpucheck::cuda_build() ? "cuda" : "cpu"; });
  m.def("device_count", &gpucheck::device_count);

  m.def(
      "check_polygons",
      [](const nb::ndarray<const i32, nb::shape<-1, 2>, nb::c_contig, nb::device::cpu>& xy,
         const In<i64, 1>& ring_start, const In<i64, 1>& shape_start, std::optional<In<i32, 1>> roles,
         const In<i32, 1>& W, const In<i32, 1>& H, i64 G) {
        const i64 n = i64(W.shape(0)), R = i64(ring_start.shape(0)) - 1;
        if (i64(H.shape(0)) != n || i64(shape_start.shape(0)) != n + 1)
          throw std::invalid_argument("W, H and shape_start must describe the same shapes");
        check_offsets(ring_start, i64(xy.shape(0)), polycore::kMaxRingVertices, "ring_start");
        check_offsets(shape_start, R, polycore::kMaxRings, "shape_start");
        polycore::check_grid(G);
        std::vector<i32> rl(size_t(R > 0 ? R : 0), -1);
        if (roles) {
          if (i64(roles->shape(0)) != R) throw std::invalid_argument("roles must have one entry per ring");
          std::copy(roles->data(), roles->data() + R, rl.begin());
        }
        for (i64 s = 0; s < n; ++s) {
          polycore::check_dims(W(s), H(s));
          polycore::check_roles(rl.data() + shape_start(s), i32(shape_start(s + 1) - shape_start(s)));
        }
        const gpucheck::Polygons in{n, xy.data(), ring_start.data(), shape_start.data(), rl.data(), W.data(),
                                    H.data(), G};
        gpucheck::PolygonResult r;
        {
          nb::gil_scoped_release release;
          r = gpucheck::check_polygons(in);
        }
        nb::dict d;
        const size_t words = r.evenodd.size();
        d["word_start"] = to_numpy(std::move(r.word_start), {size_t(n + 1)});
        d["evenodd"] = to_numpy(std::move(r.evenodd), {words});  // uint32 rows: objdel.geometry.shape.unpack_rows
        d["defects"] = to_numpy(std::move(r.defects), {size_t(n), size_t(polycore::kNumDefects)});
        return d;
      },
      nb::arg("xy").noconvert(), nb::arg("ring_start").noconvert(), nb::arg("shape_start").noconvert(),
      nb::arg("roles").noconvert().none(), nb::arg("W").noconvert(), nb::arg("H").noconvert(), nb::arg("G"));

  m.def(
      "check_masks",
      [](const In<u8, 1>& bits, const nb::ndarray<const i32, nb::shape<-1, 2>, nb::c_contig, nb::device::cpu>& hw,
         const In<i64, 1>& start) {
        const i64 n = i64(hw.shape(0));
        if (i64(start.shape(0)) != n + 1) throw std::invalid_argument("start must have one entry per mask, plus one");
        check_offsets(start, i64(bits.shape(0)), i64(polycore::kMaxDim) * polycore::kMaxDim, "start");
        std::vector<i32> H(n), W(n);
        for (i64 k = 0; k < n; ++k) {
          H[k] = hw(k, 0), W[k] = hw(k, 1);
          if (H[k] < 0 || W[k] < 0 || H[k] > polycore::kMaxDim || W[k] > polycore::kMaxDim ||
              start(k + 1) - start(k) != i64(H[k]) * W[k])
            throw std::invalid_argument("every mask must be (H, W) with 0 <= H, W <= 4096 and H * W pixels");
        }
        const gpucheck::Masks in{n, bits.data(), start.data(), H.data(), W.data()};
        gpucheck::MaskResult r;
        {
          nb::gil_scoped_release release;
          r = gpucheck::check_masks(in);
        }
        const size_t R = r.ring_area.size(), P = r.pts.size() / 2;
        nb::dict d;
        d["quads"] = to_numpy(std::move(r.quads), {size_t(n), 5});
        d["ring_start"] = to_numpy(std::move(r.ring_start), {size_t(n + 1)});
        d["pts_start"] = to_numpy(std::move(r.pts_start), {R + 1});
        d["ring_area"] = to_numpy(std::move(r.ring_area), {R});
        d["pts"] = to_numpy(std::move(r.pts), {P, 2});
        return d;
      },
      nb::arg("bits").noconvert(), nb::arg("hw").noconvert(), nb::arg("start").noconvert());
}
