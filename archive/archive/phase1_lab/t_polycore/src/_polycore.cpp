// nanobind bindings for polycore.hpp. All kernels run with the GIL released.
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include <memory>
#include <string>
#include <vector>

#include "polycore/polycore.hpp"

namespace nb = nanobind;
using namespace polycore;

using InI32 = nb::ndarray<const int32_t, nb::c_contig, nb::device::cpu>;
using InU8 = nb::ndarray<const uint8_t, nb::ndim<2>, nb::c_contig, nb::device::cpu>;

namespace {

thread_local Workspace tls_ws;  // per-thread scratch reused across calls

template <class T>
nb::ndarray<nb::numpy, T> owned(std::vector<T>&& v, std::initializer_list<size_t> shape) {
  auto* p = new std::vector<T>(std::move(v));
  nb::capsule owner(p, [](void* q) noexcept { delete static_cast<std::vector<T>*>(q); });
  return nb::ndarray<nb::numpy, T>(p->data(), shape, owner);
}

nb::ndarray<nb::numpy, bool> owned_bool(std::vector<u8>&& v, std::initializer_list<size_t> shape) {
  auto* p = new std::vector<u8>(std::move(v));
  nb::capsule owner(p, [](void* q) noexcept { delete static_cast<std::vector<u8>*>(q); });
  return nb::ndarray<nb::numpy, bool>(reinterpret_cast<bool*>(p->data()), shape, owner);
}

// Check (xy, off) and build a view. xy: int32 (N, 2) or flat (2N,); off: int32 (R+1,).
Rings make_rings(const InI32& xy, const InI32& off) {
  if (xy.ndim() == 2) {
    if (xy.shape(1) != 2) throw nb::value_error("xy must have shape (N, 2) or (2N,)");
  } else if (xy.ndim() != 1 || xy.size() % 2 != 0) {
    throw nb::value_error("xy must have shape (N, 2) or (2N,)");
  }
  if (off.ndim() != 1 || off.size() < 1) throw nb::value_error("off must be a 1-D array of R+1 ring offsets");
  if (off.size() - 1 > size_t(std::numeric_limits<i32>::max())) throw nb::value_error("too many rings");
  const size_t npts = xy.size() / 2;
  const int32_t* o = off.data();
  if (o[0] < 0) throw nb::value_error("off[0] must be >= 0");
  for (size_t i = 1; i < off.size(); ++i)
    if (o[i] < o[i - 1]) throw nb::value_error("off must be non-decreasing");
  if (size_t(o[off.size() - 1]) > npts) throw nb::value_error("off[-1] exceeds the number of vertices");
  Rings rs;
  rs.xy = xy.data();
  rs.off = o;
  rs.nrings = i32(off.size() - 1);
  return rs;
}

void check_dims(int W, int H, int64_t G) {
  if (W < 1 || H < 1 || W > kMaxDim || H > kMaxDim) throw nb::value_error("W, H must be in [1, 4096]");
  if (G < 1 || G > kMaxGrid) throw nb::value_error("G must be in [1, 2^30]");
}

std::vector<i32> check_roles(const Rings& rs, nb::handle roles_obj, bool& has) {
  has = false;
  std::vector<i32> roles;
  if (roles_obj.is_none()) return roles;
  InI32 r = nb::cast<InI32>(roles_obj);
  if (r.ndim() != 1 || r.size() != size_t(rs.nrings)) throw nb::value_error("roles must have one entry per ring");
  roles.assign(r.data(), r.data() + r.size());
  for (i32 i = 0; i < rs.nrings; ++i) {
    const i32 p = roles[i];
    if (p == -1) continue;
    if (p < 0 || p >= rs.nrings || p == i || roles[p] != -1)
      throw nb::value_error("roles[i] must be -1 (outer) or the index of an outer ring (hole)");
  }
  has = true;
  return roles;
}

const u8* check_gt(nb::handle gt_obj, int W, int H, InU8& keep) {
  if (gt_obj.is_none()) return nullptr;
  keep = nb::cast<InU8>(gt_obj);
  if (keep.shape(0) != size_t(H) || keep.shape(1) != size_t(W)) throw nb::value_error("gt must have shape (H, W)");
  return keep.data();
}

nb::object i128_to_py(i128 v) {
  // decimal string -> Python int (exact for any value)
  bool neg = v < 0;
  unsigned __int128 u = neg ? (unsigned __int128)(-(v + 1)) + 1 : (unsigned __int128)v;
  char buf[64];
  int n = 0;
  do { buf[n++] = char('0' + int(u % 10)); u /= 10; } while (u);
  std::string s = neg ? "-" : "";
  while (n) s += buf[--n];
  return nb::steal(PyLong_FromString(s.c_str(), nullptr, 10));
}

nb::dict stats_dict(const RasterStats& st) {
  nb::dict d;
  d["area_evenodd"] = st.area_evenodd;
  d["area_nonzero"] = st.area_nonzero;
  d["px_multi"] = st.px_multi;
  d["rings_mixed_sign"] = st.rings_mixed_sign;
  d["gt_area"] = st.gt_area;
  d["inter_evenodd"] = st.inter_evenodd;
  d["inter_nonzero"] = st.inter_nonzero;
  return d;
}

nb::dict validate_dict(ValidateResult& res) {
  const size_t R = size_t(res.nrings);
  nb::dict d;
  d["ring"] = owned<int32_t>(std::move(res.counts), {R, size_t(kNumDefects)});
  d["total"] = owned<int64_t>(std::vector<int64_t>(res.total.begin(), res.total.end()), {size_t(kNumDefects)});
  d["area2"] = owned<int64_t>(std::move(res.area2), {R});
  d["n_clean"] = owned<int32_t>(std::move(res.n_clean), {R});
  d["pair_cross"] = owned<int32_t>(std::move(res.pair_cross), {R, R});
  d["pair_touch"] = owned<int32_t>(std::move(res.pair_touch), {R, R});
  return d;
}

}  // namespace

NB_MODULE(_polycore, m) {
  m.doc() = "Exact integer polygon core: rasterize (even-odd + nonzero), graded validity, signed area.";

  nb::list names;
  for (int i = 0; i < kNumDefects; ++i) names.append(kDefectNames[i]);
  m.attr("DEFECTS") = nb::tuple(names);
  m.attr("MAX_RINGS_VALIDATE") = kMaxRingsValidate;

  m.def(
      "rasterize",
      [](const InI32& xy, const InI32& off, int W, int H, int64_t G, bool winding) {
        check_dims(W, H, G);
        Rings rs = make_rings(xy, off);
        const size_t npx = size_t(W) * size_t(H), R = size_t(rs.nrings);
        std::vector<u8> eo(npx), nz(npx);
        std::vector<i32> wmin(R), wmax(R), wind(winding ? npx : 0);
        RasterStats st;
        {
          nb::gil_scoped_release rel;
          RasterOut out;
          out.evenodd = eo.data();
          out.nonzero = nz.data();
          out.ring_wmin = wmin.data();
          out.ring_wmax = wmax.data();
          out.winding = winding ? wind.data() : nullptr;
          st = polycore::rasterize(rs, W, H, G, tls_ws.raster, out);
        }
        nb::dict d = stats_dict(st);
        d["evenodd"] = owned_bool(std::move(eo), {size_t(H), size_t(W)});
        d["nonzero"] = owned_bool(std::move(nz), {size_t(H), size_t(W)});
        d["ring_wmin"] = owned<int32_t>(std::move(wmin), {R});
        d["ring_wmax"] = owned<int32_t>(std::move(wmax), {R});
        if (winding) d["winding"] = owned<int32_t>(std::move(wind), {size_t(H), size_t(W)});
        return d;
      },
      nb::arg("xy").noconvert(), nb::arg("off").noconvert(), nb::arg("W"), nb::arg("H"), nb::arg("G") = 1000,
      nb::arg("winding") = false,
      "Rasterize all rings (pixel-centre sampling, half-open crossing rule). Returns a dict with "
      "bool masks 'evenodd' and 'nonzero' (H, W), per-ring 'ring_wmin'/'ring_wmax' of the ring's "
      "own winding, optional int32 'winding', and pixel statistics.");

  m.def(
      "validate",
      [](const InI32& xy, const InI32& off, nb::handle roles_obj, int64_t G, const std::string& method) {
        if (G < 1 || G > kMaxGrid) throw nb::value_error("G must be in [1, 2^30]");
        Rings rs = make_rings(xy, off);
        if (rs.nrings > kMaxRingsValidate) throw nb::value_error("too many rings for validate()");
        bool has_roles;
        std::vector<i32> roles = check_roles(rs, roles_obj, has_roles);
        PairMethod pm;
        if (method == "sweep") pm = PairMethod::Sweep;
        else if (method == "brute") pm = PairMethod::Brute;
        else throw nb::value_error("method must be 'sweep' or 'brute'");
        ValidateResult res;
        {
          nb::gil_scoped_release rel;
          polycore::validate(rs, has_roles ? roles.data() : nullptr, G, res, tls_ws.val, pm);
        }
        return validate_dict(res);
      },
      nb::arg("xy").noconvert(), nb::arg("off").noconvert(), nb::arg("roles") = nb::none(), nb::arg("G") = 1000,
      nb::arg("method") = "sweep",
      "Graded legality. roles: None (all outer) or int32 per ring, -1 = outer, p = hole of outer p. "
      "Returns dict: 'ring' (R, n_defects) counts, 'total' (n_defects,), 'area2', 'n_clean', "
      "'pair_cross', 'pair_touch' (R, R).");

  m.def(
      "signed_area2",
      [](const InI32& xy) {
        if (xy.size() % 2 != 0) throw nb::value_error("xy must hold (x, y) pairs");
        const int32_t o[2] = {0, int32_t(xy.size() / 2)};
        Rings rs;
        rs.xy = xy.data();
        rs.off = o;
        rs.nrings = 1;
        i128 a = polycore::signed_area2(rs, 0);
        return i128_to_py(a);
      },
      nb::arg("xy").noconvert(), "Twice the signed shoelace area of one ring (exact Python int); > 0 = clockwise on screen.");

  // ---- one rollout: rasterize + IoU vs gt + validate (GIL released)
  m.def(
      "score",
      [](const InI32& xy, const InI32& off, int W, int H, nb::handle roles_obj, nb::handle gt_obj, int64_t G) {
        check_dims(W, H, G);
        Rollout ro;
        ro.rings = make_rings(xy, off);
        if (ro.rings.nrings > kMaxRingsValidate) throw nb::value_error("too many rings for validate()");
        bool has_roles;
        std::vector<i32> roles = check_roles(ro.rings, roles_obj, has_roles);
        InU8 gt_keep;
        ro.gt = check_gt(gt_obj, W, H, gt_keep);
        ro.roles = has_roles ? roles.data() : nullptr;
        ro.W = W;
        ro.H = H;
        RolloutScore sc;
        {
          nb::gil_scoped_release rel;
          score_one(ro, G, tls_ws, sc);
        }
        nb::dict d = stats_dict(sc.raster);
        d["defects"] = owned<int64_t>(std::vector<int64_t>(sc.defects.begin(), sc.defects.end()), {size_t(kNumDefects)});
        return d;
      },
      nb::arg("xy").noconvert(), nb::arg("off").noconvert(), nb::arg("W"), nb::arg("H"), nb::arg("roles") = nb::none(),
      nb::arg("gt") = nb::none(), nb::arg("G") = 1000);

  // ---- many rollouts on an internal thread pool
  m.def(
      "score_batch",
      [](nb::list xys, nb::list offs, nb::list Ws, nb::list Hs, nb::object roles_list, nb::object gts,
         int64_t G, int threads) {
        const size_t n = nb::len(xys);
        if (nb::len(offs) != n || nb::len(Ws) != n || nb::len(Hs) != n)
          throw nb::value_error("xys, offs, Ws, Hs must have equal length");
        if (!roles_list.is_none() && nb::len(roles_list) != n) throw nb::value_error("roles length mismatch");
        if (!gts.is_none() && nb::len(gts) != n) throw nb::value_error("gts length mismatch");
        std::vector<Rollout> ros(n);
        std::vector<InI32> keep_xy(n), keep_off(n);
        std::vector<InU8> keep_gt(n);
        std::vector<std::vector<i32>> keep_roles(n);
        for (size_t i = 0; i < n; ++i) {
          keep_xy[i] = nb::cast<InI32>(xys[i]);
          keep_off[i] = nb::cast<InI32>(offs[i]);
          Rollout& ro = ros[i];
          ro.W = nb::cast<int>(Ws[i]);
          ro.H = nb::cast<int>(Hs[i]);
          check_dims(ro.W, ro.H, G);
          ro.rings = make_rings(keep_xy[i], keep_off[i]);
          if (ro.rings.nrings > kMaxRingsValidate) throw nb::value_error("too many rings for validate()");
          bool has_roles = false;
          if (!roles_list.is_none()) keep_roles[i] = check_roles(ro.rings, roles_list[i], has_roles);
          ro.roles = has_roles ? keep_roles[i].data() : nullptr;
          if (!gts.is_none()) ro.gt = check_gt(gts[i], ro.W, ro.H, keep_gt[i]);
        }
        std::vector<RolloutScore> out(n);
        {
          nb::gil_scoped_release rel;
          score_batch(ros.data(), i64(n), G, threads, out.data());
        }
        std::vector<int64_t> area_eo(n), area_nz(n), multi(n), gt_area(n), ieo(n), inz(n), defects(n * kNumDefects);
        std::vector<int32_t> mixed(n);
        for (size_t i = 0; i < n; ++i) {
          const RasterStats& s = out[i].raster;
          area_eo[i] = s.area_evenodd; area_nz[i] = s.area_nonzero; multi[i] = s.px_multi;
          mixed[i] = s.rings_mixed_sign; gt_area[i] = s.gt_area; ieo[i] = s.inter_evenodd; inz[i] = s.inter_nonzero;
          for (int k = 0; k < kNumDefects; ++k) defects[i * kNumDefects + k] = out[i].defects[k];
        }
        nb::dict d;
        d["area_evenodd"] = owned<int64_t>(std::move(area_eo), {n});
        d["area_nonzero"] = owned<int64_t>(std::move(area_nz), {n});
        d["px_multi"] = owned<int64_t>(std::move(multi), {n});
        d["rings_mixed_sign"] = owned<int32_t>(std::move(mixed), {n});
        d["gt_area"] = owned<int64_t>(std::move(gt_area), {n});
        d["inter_evenodd"] = owned<int64_t>(std::move(ieo), {n});
        d["inter_nonzero"] = owned<int64_t>(std::move(inz), {n});
        d["defects"] = owned<int64_t>(std::move(defects), {n, size_t(kNumDefects)});
        return d;
      },
      nb::arg("xys"), nb::arg("offs"), nb::arg("Ws"), nb::arg("Hs"), nb::arg("roles") = nb::none(),
      nb::arg("gts") = nb::none(), nb::arg("G") = 1000, nb::arg("threads") = 0,
      "Score many rollouts (rasterize + IoU counts + validate) on an internal pool of `threads` "
      "threads (0 or 1 = inline). Results are independent of the thread count.");

  m.def("hardware_threads", []() { return int(std::thread::hardware_concurrency()); });
}
