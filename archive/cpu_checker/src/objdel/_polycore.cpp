// nanobind bindings of cpp/polycore/polycore.hpp as objdel._polycore. Every kernel runs with
// the GIL released [1]; inputs are checked (caps included) before release, so bad input
// raises ValueError and never reaches the kernels. The Python API is objdel.checker.polycore.
//
// References:
//  [1] https://nanobind.readthedocs.io/en/latest/functions.html#gil-management
//  [2] https://nanobind.readthedocs.io/en/latest/ndarray.html (ndarray in/out, capsule owners)
//  [3] https://nanobind.readthedocs.io/en/latest/exceptions.html (std::invalid_argument -> ValueError)
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/vector.h>

#include <memory>
#include <vector>

#include "polycore/polycore.hpp"

namespace nb = nanobind;
using namespace polycore;

using Points = nb::ndarray<const int32_t, nb::shape<-1, 2>, nb::c_contig, nb::device::cpu>;
using Ints = nb::ndarray<const int32_t, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using Mask = nb::ndarray<const uint8_t, nb::ndim<2>, nb::c_contig, nb::device::cpu>;

namespace {

thread_local Workspace tls_ws;  // per-thread scratch reused across calls

template <class T>
nb::ndarray<nb::numpy, T> to_numpy(std::vector<T>&& v, std::initializer_list<size_t> shape) {
  auto* p = new std::vector<T>(std::move(v));
  nb::capsule owner(p, [](void* q) noexcept { delete static_cast<std::vector<T>*>(q); });
  return nb::ndarray<nb::numpy, T>(p->data(), shape, owner);
}

nb::ndarray<nb::numpy, bool> to_numpy_mask(std::unique_ptr<u8[]> m, size_t H, size_t W) {
  u8* p = m.release();
  nb::capsule owner(p, [](void* q) noexcept { delete[] static_cast<u8*>(q); });
  return nb::ndarray<nb::numpy, bool>(reinterpret_cast<bool*>(p), {H, W}, owner);
}

Rings make_rings(const Points& xy, const Ints& off) {
  if (off.shape(0) < 1) throw std::invalid_argument("off must hold R+1 ring offsets");
  Rings rs{xy.data(), off.data(), i32(off.shape(0) - 1)};
  check_rings(rs, i64(xy.shape(0)));
  return rs;
}

std::vector<i32> make_roles(const Rings& rs, const std::optional<Ints>& roles) {
  if (!roles) return {};
  if (roles->shape(0) != size_t(rs.nrings)) throw std::invalid_argument("roles must have one entry per ring");
  std::vector<i32> v(roles->data(), roles->data() + rs.nrings);
  check_roles(v.data(), rs.nrings);
  return v;
}

const u8* gt_ptr(const std::optional<Mask>& gt, i32 W, i32 H) {
  if (!gt) return nullptr;
  if (gt->shape(0) != size_t(H) || gt->shape(1) != size_t(W)) throw std::invalid_argument("gt must have shape (H, W)");
  return gt->data();
}

}  // namespace

NB_MODULE(_polycore, m) {
  m.doc() = "Exact integer polygon core (cpp/polycore/polycore.hpp); use objdel.checker.polycore.";
  nb::list names;
  for (int i = 0; i < kNumDefects; ++i) names.append(kDefectNames[i]);
  m.attr("DEFECTS") = nb::tuple(names);
  m.attr("MAX_RINGS") = kMaxRings;
  m.attr("MAX_RING_VERTICES") = kMaxRingVertices;
  m.attr("MAX_DIM") = kMaxDim;

  m.def(
      "rasterize",
      [](const Points& xy, const Ints& off, int W, int H, int64_t G) {
        check_dims(W, H);
        check_grid(G);
        const Rings rs = make_rings(xy, off);
        const size_t npx = size_t(W) * size_t(H), R = size_t(rs.nrings);
        std::unique_ptr<u8[]> eo(new u8[npx]), nz(new u8[npx]);  // zeroed by rasterize
        std::vector<i32> wmin(R), wmax(R);
        RasterStats st;
        {
          nb::gil_scoped_release release;
          st = rasterize(rs, W, H, G, tls_ws.raster, RasterOut{eo.get(), nz.get(), wmin.data(), wmax.data()});
        }
        nb::dict d;
        d["evenodd"] = to_numpy_mask(std::move(eo), H, W);
        d["nonzero"] = to_numpy_mask(std::move(nz), H, W);
        d["ring_wmin"] = to_numpy(std::move(wmin), {R});
        d["ring_wmax"] = to_numpy(std::move(wmax), {R});
        d["area_evenodd"] = st.area_evenodd;
        d["area_nonzero"] = st.area_nonzero;
        d["px_multi"] = st.px_multi;
        d["rings_mixed_sign"] = st.rings_mixed_sign;
        return d;
      },
      nb::arg("xy").noconvert(), nb::arg("off").noconvert(), nb::arg("W"), nb::arg("H"), nb::arg("G"));

  m.def(
      "validate",
      [](const Points& xy, const Ints& off, std::optional<Ints> roles, int64_t G) {
        check_grid(G);
        const Rings rs = make_rings(xy, off);
        const std::vector<i32> rl = make_roles(rs, roles);
        ValidateResult res;
        {
          nb::gil_scoped_release release;
          validate(rs, rl.empty() ? nullptr : rl.data(), G, res, tls_ws.val);
        }
        const size_t R = size_t(res.nrings);
        nb::dict d;
        d["total"] = to_numpy(std::vector<int64_t>(res.total.begin(), res.total.end()), {size_t(kNumDefects)});
        d["per_ring"] = to_numpy(std::move(res.counts), {R, size_t(kNumDefects)});
        d["area2"] = to_numpy(std::move(res.area2), {R});
        d["n_clean"] = to_numpy(std::move(res.n_clean), {R});
        d["pair_cross"] = to_numpy(std::move(res.pair_cross), {R, R});
        d["pair_touch"] = to_numpy(std::move(res.pair_touch), {R, R});
        return d;
      },
      nb::arg("xy").noconvert(), nb::arg("off").noconvert(), nb::arg("roles").noconvert().none(), nb::arg("G"));

  m.def(
      "signed_area2",
      [](const Points& xy) {
        const i32 off[2] = {0, i32(xy.shape(0))};
        const i128 a = signed_area2(Rings{xy.data(), off, 1}, 0);
        return nb::make_tuple(int64_t(a >> 64), uint64_t(a));  // exact: a = hi * 2^64 + lo
      },
      nb::arg("xy").noconvert());

  m.def(
      "score_batch",
      [](const std::vector<Points>& xys, const std::vector<Ints>& offs, const std::vector<int>& Ws,
         const std::vector<int>& Hs, std::optional<std::vector<std::optional<Ints>>> roles,
         std::optional<std::vector<std::optional<Mask>>> gts, int64_t G, int threads) {
        const size_t n = xys.size();
        if (offs.size() != n || Ws.size() != n || Hs.size() != n || (roles && roles->size() != n) ||
            (gts && gts->size() != n))
          throw std::invalid_argument("score_batch: all per-rollout lists must have the same length");
        check_grid(G);
        std::vector<Rollout> ros(n);
        std::vector<std::vector<i32>> keep_roles(n);
        for (size_t i = 0; i < n; ++i) {
          check_dims(Ws[i], Hs[i]);
          ros[i].rings = make_rings(xys[i], offs[i]);
          ros[i].W = Ws[i];
          ros[i].H = Hs[i];
          if (roles) keep_roles[i] = make_roles(ros[i].rings, (*roles)[i]);
          ros[i].roles = keep_roles[i].empty() ? nullptr : keep_roles[i].data();
          if (gts) ros[i].gt = gt_ptr((*gts)[i], Ws[i], Hs[i]);
        }
        std::vector<RolloutScore> out(n);
        {
          nb::gil_scoped_release release;
          score_batch(ros.data(), i64(n), G, threads, out.data());
        }
        std::vector<int64_t> aeo(n), anz(n), multi(n), gta(n), ieo(n), inz(n), defects(n * kNumDefects);
        std::vector<int32_t> mixed(n);
        for (size_t i = 0; i < n; ++i) {
          const RasterStats& s = out[i].raster;
          aeo[i] = s.area_evenodd, anz[i] = s.area_nonzero, multi[i] = s.px_multi, mixed[i] = s.rings_mixed_sign;
          gta[i] = s.gt_area, ieo[i] = s.inter_evenodd, inz[i] = s.inter_nonzero;
          std::copy(out[i].defects.begin(), out[i].defects.end(), defects.begin() + i * kNumDefects);
        }
        nb::dict d;
        d["area_evenodd"] = to_numpy(std::move(aeo), {n});
        d["area_nonzero"] = to_numpy(std::move(anz), {n});
        d["px_multi"] = to_numpy(std::move(multi), {n});
        d["rings_mixed_sign"] = to_numpy(std::move(mixed), {n});
        d["gt_area"] = to_numpy(std::move(gta), {n});
        d["inter_evenodd"] = to_numpy(std::move(ieo), {n});
        d["inter_nonzero"] = to_numpy(std::move(inz), {n});
        d["defects"] = to_numpy(std::move(defects), {n, size_t(kNumDefects)});
        return d;
      },
      nb::arg("xys"), nb::arg("offs"), nb::arg("Ws"), nb::arg("Hs"), nb::arg("roles").none(),
      nb::arg("gts").none(), nb::arg("G"), nb::arg("threads"));
}
