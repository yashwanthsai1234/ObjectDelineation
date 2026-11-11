// nanobind binding of cuda/gpucheck.h as objdel._gpucheck (CPU build by default; `make cuda` builds it with nvcc).
// Inputs are the packed arrays of data/checker_eval (experiments/cuda_checker/prepare_data.py); they are checked here
// with gpucheck.h's caps before the GIL is released, so bad input raises ValueError and never reaches the kernels.
// The GIL is released during the checker call (other Python threads run; gpucheck's mutex orders concurrent calls)
// and taken back before the result arrays are made. References: https://nanobind.readthedocs.io/en/latest/ndarray.html
// , https://nanobind.readthedocs.io/en/latest/functions.html#gil-management .
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>

#include <cmath>
#include <optional>
#include <stdexcept>
#include <vector>

#include "gpucheck.h"

namespace nb = nanobind;
using gpucheck::i32;
using gpucheck::i64;
using gpucheck::u8;

template <class T, int D>
using In = nb::ndarray<const T, nb::ndim<D>, nb::c_contig, nb::device::cpu>;
using XY = nb::ndarray<const i32, nb::shape<-1, 2>, nb::c_contig, nb::device::cpu>;

namespace {

template <class T>
nb::ndarray<nb::numpy, T> to_numpy(std::vector<T>&& v, std::initializer_list<size_t> shape) {
  auto* p = new std::vector<T>(std::move(v));
  nb::capsule owner(p, [](void* q) noexcept { delete static_cast<std::vector<T>*>(q); });
  return nb::ndarray<nb::numpy, T>(p->data(), shape, owner);
}

nb::ndarray<nb::numpy, bool> to_bool(std::vector<u8>&& v) {  // 0 / 1 bytes as a numpy bool array
  auto* p = new std::vector<u8>(std::move(v));
  nb::capsule owner(p, [](void* q) noexcept { delete static_cast<std::vector<u8>*>(q); });
  return nb::ndarray<nb::numpy, bool>(reinterpret_cast<bool*>(p->data()), {p->size()}, owner);
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

// Checked polygon input (roles defaulted to all outers, kept alive here).
struct PolygonsIn {
  std::vector<i32> roles;
  gpucheck::Polygons in;
  PolygonsIn(const XY& xy, const In<i64, 1>& ring_start, const In<i64, 1>& shape_start,
             const std::optional<In<i32, 1>>& rl, const In<i32, 1>& W, const In<i32, 1>& H, i64 G) {
    const i64 n = i64(W.shape(0)), R = i64(ring_start.shape(0)) - 1;
    if (i64(H.shape(0)) != n || i64(shape_start.shape(0)) != n + 1)
      throw std::invalid_argument("W, H and shape_start must describe the same shapes");
    check_offsets(ring_start, i64(xy.shape(0)), gpucheck::kMaxRingVertices, "ring_start");
    check_offsets(shape_start, R, gpucheck::kMaxRings, "shape_start");
    gpucheck::check_grid(G);
    roles.assign(size_t(R > 0 ? R : 0), -1);
    if (rl) {
      if (i64(rl->shape(0)) != R) throw std::invalid_argument("roles must have one entry per ring");
      std::copy(rl->data(), rl->data() + R, roles.begin());
    }
    for (i64 s = 0; s < n; ++s) {
      gpucheck::check_dims(W(s), H(s));
      gpucheck::check_roles(roles.data() + shape_start(s), i32(shape_start(s + 1) - shape_start(s)));
    }
    in = {n, xy.data(), ring_start.data(), shape_start.data(), roles.data(), W.data(), H.data(), G};
  }
};

// Checked masks: (H, W) per mask, H * W bytes each.
struct MasksIn {
  std::vector<i32> H, W;
  gpucheck::Masks in;
  MasksIn(const In<u8, 1>& bits, const XY& hw, const In<i64, 1>& start) {
    const i64 n = i64(hw.shape(0));
    if (i64(start.shape(0)) != n + 1) throw std::invalid_argument("start must have one entry per mask, plus one");
    check_offsets(start, i64(bits.shape(0)), i64(gpucheck::kMaxDim) * gpucheck::kMaxDim, "start");
    H.resize(size_t(n)), W.resize(size_t(n));
    for (i64 k = 0; k < n; ++k) {
      H[size_t(k)] = hw(k, 0), W[size_t(k)] = hw(k, 1);
      if (H[size_t(k)] < 0 || W[size_t(k)] < 0 || H[size_t(k)] > gpucheck::kMaxDim ||
          W[size_t(k)] > gpucheck::kMaxDim || start(k + 1) - start(k) != i64(H[size_t(k)]) * W[size_t(k)])
        throw std::invalid_argument("every mask must be (H, W) with 0 <= H, W <= 4096 and H * W pixels");
    }
    in = {n, bits.data(), start.data(), H.data(), W.data()};
  }
};

nb::dict scores_dict(gpucheck::Scores&& r, i64 n) {
  nb::dict d;
  const size_t S = r.sig.size();
  d["defects"] = to_numpy(std::move(r.defects), {size_t(n), size_t(gpucheck::kNumDefects)});
  d["area"] = to_numpy(std::move(r.area), {size_t(n)});
  d["sum_x"] = to_numpy(std::move(r.sum_x), {size_t(n)});
  d["sum_y"] = to_numpy(std::move(r.sum_y), {size_t(n)});
  d["iou_tc"] = to_numpy(std::move(r.iou_tc), {size_t(n)});
  d["f5"] = to_numpy(std::move(r.f5), {size_t(n)});
  d["g_area"] = to_numpy(std::move(r.g_area), {size_t(n)});
  d["topo_freq"] = to_numpy(std::move(r.topo_freq), {size_t(n)});
  d["sig_start"] = to_numpy(std::move(r.sig_start), {size_t(n + 1)});
  d["sig"] = to_numpy(std::move(r.sig), {S});
  return d;
}

void check_tables(const In<i64, 1>& tables, i64 n) {
  if (i64(tables.shape(0)) != n) throw std::invalid_argument("tables must have one handle per object");
}

// add_tables' arguments, checked (frame and rater caps, shapes, the dist precondition, the scalars).
gpucheck::Tables tables_of(const In<u8, 3>& raters, const In<u8, 2>& uni, const In<u8, 2>& core, const In<u8, 2>& vote,
                           const In<float, 2>& dist, const In<i64, 1>& rater_area, const In<u8, 1>& voter,
                           const In<i64, 1>& sig_start, const In<i32, 1>& sig, i64 n_vote, i64 perimeter, double t5,
                           i32 min_area) {
  const i64 k = i64(raters.shape(0)), H = i64(raters.shape(1)), W = i64(raters.shape(2));
  gpucheck::check_dims(W, H);
  if (k < 2 || k > 32) throw std::invalid_argument("tables need 2 to 32 raters");
  for (const auto* a : {&uni, &core, &vote})
    if (i64(a->shape(0)) != H || i64(a->shape(1)) != W) throw std::invalid_argument("planes must be (H, W)");
  if (i64(dist.shape(0)) != H || i64(dist.shape(1)) != W) throw std::invalid_argument("dist must be (H, W)");
  if (i64(rater_area.shape(0)) != k || i64(voter.shape(0)) != k || i64(sig_start.shape(0)) != k + 1)
    throw std::invalid_argument("rater_area, voter and sig_start must have one entry per rater (+1)");
  check_offsets(sig_start, i64(sig.shape(0)), i64(H) * W, "sig_start");
  for (i64 i = 0; i < H * W; ++i) {  // the exactness of the band-depth sums (gpucheck.cpp, scoring)
    const float d = dist.data()[i];
    if (!(d == 0.0f || (d >= 1.0f && d <= 16777216.0f)))
      throw std::invalid_argument("dist values must be 0 or in [1, 2^24]");
  }
  if (perimeter < 1 || !(t5 > 0) || !std::isfinite(t5) || min_area < 0)
    throw std::invalid_argument("perimeter >= 1, t5 > 0 and min_area >= 0 are required");
  return {i32(H),           i32(W),       i32(k),           raters.data(), uni.data(), core.data(),
          vote.data(),      dist.data(),  rater_area.data(), voter.data(), sig_start.data(), sig.data(),
          n_vote,           perimeter,    t5,               min_area};
}

struct TablesIn {  // add_tables_many: one item's arrays, kept alive during the call
  In<u8, 3> raters;
  In<u8, 2> uni, core, vote;
  In<float, 2> dist;
  In<i64, 1> rater_area;
  In<u8, 1> voter;
  In<i64, 1> sig_start;
  In<i32, 1> sig;
};

// timing=True: the call's per-stage milliseconds (stage -> ms, in run order) as the result's "timing".
nb::dict times_dict(const gpucheck::Times& times) {
  nb::dict t;
  for (const auto& [name, ms] : times) t[name.c_str()] = ms;
  return t;
}

}  // namespace

NB_MODULE(_gpucheck, m) {
  m.doc() = "Batched exact checker (cuda/gpucheck.h): polygon legality (DEFECTS) and even-odd masks "
            "(bit-packed rows), polygonize's rings, and the reward's scorer against expert tables resident on the device.";
  m.def("build", [] { return gpucheck::cuda_build() ? "cuda" : "cpu"; });
  m.def("device_count", &gpucheck::device_count);
  m.def("stats", [] {  // gpucheck.h Stats
    const gpucheck::Stats s = gpucheck::stats();
    nb::dict d;
    d["waits"] = s.waits, d["drain_syncs"] = s.drain_syncs, d["label_pairs"] = s.label_pairs;
    d["sync_alloc"] = s.sync_alloc, d["sync_free"] = s.sync_free;
    d["pool_allocs"] = s.pool_allocs, d["pool_frees"] = s.pool_frees, d["pinned_grows"] = s.pinned_grows;
    d["arena_bytes"] = s.arena_bytes, d["table_bytes"] = s.table_bytes;
    d["live_device_bytes"] = s.live_device_bytes;
    d["pool_reserved"] = s.pool_reserved, d["default_pool_reserved"] = s.default_pool_reserved;
    d["eval_scratch_peak"] = s.eval_scratch_peak;
    return d;
  });
  m.def("reset_stats", &gpucheck::reset_stats);
  m.def("trim", &gpucheck::trim, nb::arg("bytes"));
  nb::list names;
  for (const char* name : gpucheck::kDefectNames) names.append(name);
  m.attr("DEFECTS") = nb::tuple(names);

  m.def(
      "check_polygons",
      [](const XY& xy, const In<i64, 1>& ring_start, const In<i64, 1>& shape_start, std::optional<In<i32, 1>> roles,
         const In<i32, 1>& W, const In<i32, 1>& H, i64 G, bool timing) {
        const PolygonsIn p(xy, ring_start, shape_start, roles, W, H, G);
        const i64 n = p.in.n;
        gpucheck::PolygonResult r;
        gpucheck::Times times;
        {
          nb::gil_scoped_release release;
          r = gpucheck::check_polygons(p.in, timing ? &times : nullptr);
        }
        nb::dict d;
        const size_t words = r.evenodd.size();
        d["word_start"] = to_numpy(std::move(r.word_start), {size_t(n + 1)});
        d["evenodd"] = to_numpy(std::move(r.evenodd), {words});  // uint32 rows: objdel.geometry.shape.unpack_rows
        d["defects"] = to_numpy(std::move(r.defects), {size_t(n), size_t(gpucheck::kNumDefects)});
        if (timing) d["timing"] = times_dict(times);
        return d;
      },
      nb::arg("xy").noconvert(), nb::arg("ring_start").noconvert(), nb::arg("shape_start").noconvert(),
      nb::arg("roles").noconvert().none(), nb::arg("W").noconvert(), nb::arg("H").noconvert(), nb::arg("G"),
      nb::kw_only(), nb::arg("timing") = false);

  m.def("test_hooks", [] {  // () in every shipped build; tests/variants.py checks its builds with it
    nb::list h;
    for (const std::string& s : gpucheck::test_hooks()) h.append(s);
    return nb::tuple(h);
  });

  // chunk_pixels, max_scratch_bytes and seed_rows are for tests (gpucheck.h MaskOptions; defaults = production).
  m.def(
      "check_masks",
      [](const In<u8, 1>& bits, const XY& hw, const In<i64, 1>& start, bool timing, i64 chunk_pixels,
         i64 max_scratch_bytes, i32 seed_rows) {
        const MasksIn k(bits, hw, start);
        const i64 n = k.in.n;
        gpucheck::MaskResult r;
        gpucheck::Times times;
        {
          nb::gil_scoped_release release;
          r = gpucheck::check_masks(k.in, timing ? &times : nullptr, {seed_rows, chunk_pixels, max_scratch_bytes});
        }
        const size_t R = r.ring_area.size(), P = r.pts.size() / 2;
        nb::dict d;
        d["quads"] = to_numpy(std::move(r.quads), {size_t(n), 5});
        d["ring_start"] = to_numpy(std::move(r.ring_start), {size_t(n + 1)});
        d["pts_start"] = to_numpy(std::move(r.pts_start), {R + 1});
        d["ring_area"] = to_numpy(std::move(r.ring_area), {R});
        d["pts"] = to_numpy(std::move(r.pts), {P, 2});
        d["ring_part"] = to_numpy(std::move(r.ring_part), {R});
        d["betti"] = to_numpy(std::move(r.betti), {size_t(n), 2});
        if (timing) d["timing"] = times_dict(times);
        return d;
      },
      nb::arg("bits").noconvert(), nb::arg("hw").noconvert(), nb::arg("start").noconvert(), nb::kw_only(),
      nb::arg("timing") = false, nb::arg("chunk_pixels") = gpucheck::MaskOptions{}.chunk_pixels,
      nb::arg("max_scratch_bytes") = gpucheck::MaskOptions{}.max_scratch_bytes,
      nb::arg("seed_rows") = gpucheck::MaskOptions{}.seed_rows);

  // One ExpertCache's tables to the device -> handle (objdel.reward.ExpertCache.device).
  m.def(
      "add_tables",
      [](const In<u8, 3>& raters, const In<u8, 2>& uni, const In<u8, 2>& core, const In<u8, 2>& vote,
         const In<float, 2>& dist, const In<i64, 1>& rater_area, const In<u8, 1>& voter, const In<i64, 1>& sig_start,
         const In<i32, 1>& sig, i64 n_vote, i64 perimeter, double t5, i32 min_area) {
        const gpucheck::Tables t = tables_of(raters, uni, core, vote, dist, rater_area, voter, sig_start, sig, n_vote,
                                             perimeter, t5, min_area);
        nb::gil_scoped_release release;
        return gpucheck::add_tables(t);
      },
      nb::arg("raters").noconvert(), nb::arg("uni").noconvert(), nb::arg("core").noconvert(),
      nb::arg("vote").noconvert(), nb::arg("dist").noconvert(), nb::arg("rater_area").noconvert(),
      nb::arg("voter").noconvert(), nb::arg("sig_start").noconvert(), nb::arg("sig").noconvert(), nb::arg("n_vote"),
      nb::arg("perimeter"), nb::arg("t5"), nb::arg("min_area"));
  // Many caches' tables in one batch (one host wait) -> int64 handles; each item is add_tables' 13 arguments. Every
  // item is checked before anything is uploaded.
  m.def(
      "add_tables_many",
      [](const nb::list& items) {
        std::vector<TablesIn> keep;  // the arrays stay referenced during the call
        std::vector<gpucheck::Tables> ts;
        for (nb::handle item : items) {
          const nb::tuple a = nb::cast<nb::tuple>(item);
          if (a.size() != 13) throw std::invalid_argument("each item must be the 13 arguments of add_tables");
          keep.push_back(TablesIn{nb::cast<In<u8, 3>>(a[0], false),  nb::cast<In<u8, 2>>(a[1], false),
                                  nb::cast<In<u8, 2>>(a[2], false),  nb::cast<In<u8, 2>>(a[3], false),
                                  nb::cast<In<float, 2>>(a[4], false), nb::cast<In<i64, 1>>(a[5], false),
                                  nb::cast<In<u8, 1>>(a[6], false),  nb::cast<In<i64, 1>>(a[7], false),
                                  nb::cast<In<i32, 1>>(a[8], false)});
          const TablesIn& k = keep.back();
          ts.push_back(tables_of(k.raters, k.uni, k.core, k.vote, k.dist, k.rater_area, k.voter, k.sig_start, k.sig,
                                 nb::cast<i64>(a[9]), nb::cast<i64>(a[10]), nb::cast<double>(a[11]),
                                 nb::cast<i32>(a[12])));
        }
        std::vector<i64> h;
        {
          nb::gil_scoped_release release;
          h = gpucheck::add_tables_many(ts);
        }
        const size_t m = h.size();
        return to_numpy(std::move(h), {m});
      },
      nb::arg("items"));
  // Many caches' tables built on the device from their rater planes (objdel.reward.cache.device_many): one upload and
  // one host wait per batch. Every item is checked before anything is uploaded.
  m.def(
      "add_tables_from_raters",
      [](const nb::list& raters, const In<double, 1>& t5, const In<i64, 1>& min_area) {
        const i64 m = i64(raters.size());
        if (i64(t5.shape(0)) != m || i64(min_area.shape(0)) != m)
          throw std::invalid_argument("t5 and min_area must have one entry per cache");
        std::vector<In<u8, 3>> keep;  // the arrays stay referenced during the call
        std::vector<gpucheck::Raters> in;
        for (i64 c = 0; c < m; ++c) {
          keep.push_back(nb::cast<In<u8, 3>>(raters[size_t(c)], false));
          const In<u8, 3>& a = keep.back();
          const i64 k = i64(a.shape(0)), H = i64(a.shape(1)), W = i64(a.shape(2));
          gpucheck::check_dims(W, H);
          if (k < 2 || k > 32) throw std::invalid_argument("tables need 2 to 32 raters");
          if (!(t5(c) > 0) || !std::isfinite(t5(c)) || min_area(c) < 0 || min_area(c) > INT32_MAX)
            throw std::invalid_argument("t5 > 0 and 0 <= min_area < 2^31 are required");
          in.push_back(gpucheck::Raters{i32(H), i32(W), i32(k), a.data(), t5(c), i32(min_area(c))});
        }
        gpucheck::BuiltTables r;
        {
          nb::gil_scoped_release release;
          r = gpucheck::add_tables_from_raters(in);
        }
        const size_t K = r.voter.size(), S = r.sig.size();
        nb::dict d;
        d["handles"] = to_numpy(std::move(r.handles), {size_t(m)});
        d["perimeter"] = to_numpy(std::move(r.perimeter), {size_t(m)});
        d["rater_area_start"] = to_numpy(std::move(r.rater_area_start), {size_t(m + 1)});
        d["rater_area"] = to_numpy(std::move(r.rater_area), {K});
        d["sig_start"] = to_numpy(std::move(r.sig_start), {K + 1});
        d["sig"] = to_numpy(std::move(r.sig), {S});
        d["voter"] = to_numpy(std::move(r.voter), {K});
        d["n_vote"] = to_numpy(std::move(r.n_vote), {size_t(m)});
        return d;
      },
      nb::arg("raters"), nb::arg("t5").noconvert(), nb::arg("min_area").noconvert());
  // objdel.reward.cache.edt: (H, W) uint8, nonzero = not a source -> float32 (H, W).
  m.def(
      "edt",
      [](const In<u8, 2>& nonzero) {
        const i64 H = i64(nonzero.shape(0)), W = i64(nonzero.shape(1));
        if (H > gpucheck::kMaxEdtDim || W > gpucheck::kMaxEdtDim)
          throw std::invalid_argument("edt: the map must be (H, W) with 0 <= H, W <= 32767");
        std::vector<float> r;
        {
          nb::gil_scoped_release release;
          r = gpucheck::edt(nonzero.data(), i32(H), i32(W));
        }
        return to_numpy(std::move(r), {size_t(H), size_t(W)});
      },
      nb::arg("nonzero").noconvert());
  m.def("drop_tables", &gpucheck::drop_tables, nb::arg("handle"));
  m.def("resident", &gpucheck::resident);

  m.def(
      "score_polygons",
      [](const XY& xy, const In<i64, 1>& ring_start, const In<i64, 1>& shape_start, std::optional<In<i32, 1>> roles,
         const In<i32, 1>& W, const In<i32, 1>& H, i64 G, const In<i64, 1>& tables, bool timing) {
        const PolygonsIn p(xy, ring_start, shape_start, roles, W, H, G);
        check_tables(tables, p.in.n);
        gpucheck::Scores r;
        gpucheck::Times times;
        {
          nb::gil_scoped_release release;
          r = gpucheck::score_polygons(p.in, tables.data(), timing ? &times : nullptr);
        }
        nb::dict d = scores_dict(std::move(r), p.in.n);
        if (timing) d["timing"] = times_dict(times);
        return d;
      },
      nb::arg("xy").noconvert(), nb::arg("ring_start").noconvert(), nb::arg("shape_start").noconvert(),
      nb::arg("roles").noconvert().none(), nb::arg("W").noconvert(), nb::arg("H").noconvert(), nb::arg("G"),
      nb::arg("tables").noconvert(), nb::kw_only(), nb::arg("timing") = false);

  m.def(
      "score_masks",
      [](const In<u8, 1>& bits, const XY& hw, const In<i64, 1>& start, const In<i64, 1>& tables, bool timing) {
        const MasksIn k(bits, hw, start);
        check_tables(tables, k.in.n);
        gpucheck::Scores r;
        gpucheck::Times times;
        {
          nb::gil_scoped_release release;
          r = gpucheck::score_masks(k.in, tables.data(), timing ? &times : nullptr);
        }
        nb::dict d = scores_dict(std::move(r), k.in.n);
        if (timing) d["timing"] = times_dict(times);
        return d;
      },
      nb::arg("bits").noconvert(), nb::arg("hw").noconvert(), nb::arg("start").noconvert(),
      nb::arg("tables").noconvert(), nb::kw_only(), nb::arg("timing") = false);

  // objdel.checker.topology on the engine (gpucheck.h Topology): per-component numbers; pixels only with want_masks.
  m.def(
      "topology",
      [](const In<u8, 1>& bits, const XY& hw, const In<i64, 1>& start, const In<i64, 1>& min_area, bool want_masks) {
        const MasksIn k(bits, hw, start);
        const i64 n = k.in.n;
        if (i64(min_area.shape(0)) != n) throw std::invalid_argument("min_area must have one entry per mask");
        for (i64 i = 0; i < n; ++i)
          if (min_area(i) < 0) throw std::invalid_argument("min_area must be >= 0");
        gpucheck::Topology r;
        {
          nb::gil_scoped_release release;
          r = gpucheck::topology(k.in, min_area.data(), want_masks);
        }
        const size_t C = r.comp_area.size(), Hn = r.hole_area.size(), S = r.sig.size(), N = r.masks.size();
        nb::dict d;
        d["comp_start"] = to_numpy(std::move(r.comp_start), {size_t(n + 1)});
        d["comp_area"] = to_numpy(std::move(r.comp_area), {C});
        d["comp_holes"] = to_numpy(std::move(r.comp_holes), {C});
        d["hole_start"] = to_numpy(std::move(r.hole_start), {C + 1});
        d["hole_area"] = to_numpy(std::move(r.hole_area), {Hn});
        d["sig_start"] = to_numpy(std::move(r.sig_start), {size_t(n + 1)});
        d["sig"] = to_numpy(std::move(r.sig), {S});
        if (want_masks) d["masks"] = to_numpy(std::move(r.masks), {N});
        return d;
      },
      nb::arg("bits").noconvert(), nb::arg("hw").noconvert(), nb::arg("start").noconvert(),
      nb::arg("min_area").noconvert(), nb::arg("want_masks") = false);

  // objdel.eval.device: the item-level counts of metrics() for sets of K sample polygons (gpucheck.h EvalSets).
  m.def(
      "eval_sets",
      [](const XY& xy, const In<i64, 1>& ring_start, const In<i64, 1>& shape_start, std::optional<In<i32, 1>> roles,
         const In<i32, 1>& W, const In<i32, 1>& H, i64 G, const In<i64, 1>& tables, i32 K, const In<u8, 1>& legal,
         const In<i32, 3>& soft_cuts, const In<i64, 1>& nsd_cut, const In<i32, 1>& band_d, i64 max_scratch_bytes) {
        const PolygonsIn p(xy, ring_start, shape_start, roles, W, H, G);
        const i64 n = i64(tables.shape(0)), T = i64(soft_cuts.shape(1));
        if (K < 1 || K > 256) throw std::invalid_argument("eval_sets: K must be in [1, 256]");
        if (p.in.n != n * K || i64(legal.shape(0)) != n * K)
          throw std::invalid_argument("eval_sets: K shapes and K legal flags per set");
        if (i64(soft_cuts.shape(0)) != n || i64(soft_cuts.shape(2)) != 2 || T < 1 || T > gpucheck::kMaxThresholds)
          throw std::invalid_argument("eval_sets: soft_cuts must be (sets, 1..16 thresholds, 2)");
        if (i64(nsd_cut.shape(0)) != n || i64(band_d.shape(0)) != n)
          throw std::invalid_argument("eval_sets: one nsd_cut and one band_d per set");
        if (max_scratch_bytes < 1) throw std::invalid_argument("eval_sets: max_scratch_bytes must be >= 1");
        const gpucheck::EvalSets es{n,           K,           i32(T), tables.data(), legal.data(), soft_cuts.data(),
                                    nsd_cut.data(), band_d.data(), max_scratch_bytes};
        for (i64 s = 0; s < n; ++s) {
          if (nsd_cut(s) < 0 || band_d(s) < 1 || band_d(s) > 31)
            throw std::invalid_argument("eval_sets: nsd_cut >= 0 and 1 <= band_d <= 31 are required");
          for (i64 t = 0; t < T; ++t)
            if (soft_cuts(s, t, 0) < 0 || soft_cuts(s, t, 0) > K || soft_cuts(s, t, 1) < 0)
              throw std::invalid_argument("eval_sets: cuts must be in [0, K] and [0, k]");
        }
        gpucheck::EvalCounts r;
        {
          nb::gil_scoped_release release;
          r = gpucheck::eval_sets(p.in, es);
        }
        const size_t N = size_t(n), k = size_t(r.k), M = size_t(K) + k + 1, Kz = size_t(K);
        std::vector<i64> area_raw(N * Kz);
        for (size_t q = 0; q < N * Kz; ++q) area_raw[q] = r.inter_raw[(q / Kz * Kz + q % Kz) * Kz + q % Kz];
        nb::dict d;
        d["inter"] = to_numpy(std::move(r.inter), {N, M, M});
        d["inter_raw"] = to_numpy(std::move(r.inter_raw), {N, Kz, Kz});
        d["area_raw"] = to_numpy(std::move(area_raw), {N, Kz});
        d["band_inter"] = to_numpy(std::move(r.band_inter), {N, M, M});
        d["soft"] = to_numpy(std::move(r.soft), {N, size_t(T), 3});
        d["surf"] = to_numpy(std::move(r.surf), {N, Kz + k, 3});
        d["nsd_out"] = to_numpy(std::move(r.nsd), {N, Kz, k, 2, 3});
        d["vote_iou_tc"] = to_numpy(std::move(r.vote_iou_tc), {N});
        return d;
      },
      nb::arg("xy").noconvert(), nb::arg("ring_start").noconvert(), nb::arg("shape_start").noconvert(),
      nb::arg("roles").noconvert().none(), nb::arg("W").noconvert(), nb::arg("H").noconvert(), nb::arg("G"),
      nb::arg("tables").noconvert(), nb::arg("K"), nb::arg("legal").noconvert(), nb::arg("soft_cuts").noconvert(),
      nb::arg("nsd_cut").noconvert(), nb::arg("band_d").noconvert(), nb::kw_only(),
      nb::arg("max_scratch_bytes") = i64(256) << 20);

  // The target builder's ladder: every rung of one name against its reference mask in one call (gpucheck.h).
  m.def(
      "certify_polygons",
      [](const XY& xy, const In<i64, 1>& ring_start, const In<i64, 1>& shape_start, std::optional<In<i32, 1>> roles,
         const In<i32, 1>& W, const In<i32, 1>& H, i64 G, const In<u8, 1>& ref_bits, const XY& ref_hw,
         const In<i64, 1>& ref_start, const In<i32, 1>& ref_of, const In<i64, 1>& min_area,
         std::optional<XY> pairs, std::optional<XY> ref_pairs, std::optional<In<i64, 1>> pair_min_area) {
        const PolygonsIn p(xy, ring_start, shape_start, roles, W, H, G);
        const MasksIn r(ref_bits, ref_hw, ref_start);
        const i64 n = p.in.n, nr = r.in.n;
        if (i64(ref_of.shape(0)) != n || i64(min_area.shape(0)) != n)
          throw std::invalid_argument("ref_of and min_area must have one entry per shape");
        for (i64 s = 0; s < n; ++s) {
          if (ref_of(s) < 0 || ref_of(s) >= nr) throw std::invalid_argument("ref_of must name a reference plane");
          if (r.H[size_t(ref_of(s))] != H(s) || r.W[size_t(ref_of(s))] != W(s))
            throw std::invalid_argument("a shape's reference must have the shape's frame (H, W)");
          if (min_area(s) < 0) throw std::invalid_argument("min_area must be >= 0");
        }
        gpucheck::Pairs q;
        if (pairs || ref_pairs || pair_min_area) {
          if (!pairs || !ref_pairs || !pair_min_area)
            throw std::invalid_argument("pairs, ref_pairs and pair_min_area go together");
          q = {i64(pairs->shape(0)), pairs->data(), ref_pairs->data(), pair_min_area->data()};
          if (i64(ref_pairs->shape(0)) != q.n || i64(pair_min_area->shape(0)) != q.n)
            throw std::invalid_argument("pairs, ref_pairs and pair_min_area must have one entry per pair");
          for (i64 k = 0; k < q.n; ++k) {
            const i32 a = q.shapes[2 * k], b = q.shapes[2 * k + 1], ra = q.refs[2 * k], rb = q.refs[2 * k + 1];
            if (a < 0 || b < 0 || a >= n || b >= n || ra < 0 || rb < 0 || ra >= nr || rb >= nr)
              throw std::invalid_argument("pairs must name shapes and ref_pairs reference planes of this call");
            if (W(a) != W(b) || H(a) != H(b) || r.W[size_t(ra)] != r.W[size_t(rb)] || r.H[size_t(ra)] != r.H[size_t(rb)])
              throw std::invalid_argument("the two planes of a pair must share one frame");
            if (q.min_area[k] < 0) throw std::invalid_argument("pair_min_area must be >= 0");
          }
        }
        gpucheck::Certified c;
        {
          nb::gil_scoped_release release;
          c = gpucheck::certify_polygons(p.in, r.in, ref_of.data(), min_area.data(), q);
        }
        nb::dict d;
        const size_t S = c.sig.size(), S2 = c.ref_sig.size(), Q = c.pair_sig.size(), Q2 = c.pair_ref_sig.size();
        d["legal"] = to_bool(std::move(c.legal));
        d["inter"] = to_numpy(std::move(c.inter), {size_t(n)});
        d["union"] = to_numpy(std::move(c.uni), {size_t(n)});
        d["equal"] = to_bool(std::move(c.equal));
        d["sig_start"] = to_numpy(std::move(c.sig_start), {size_t(n + 1)});
        d["sig"] = to_numpy(std::move(c.sig), {S});
        d["ref_sig_start"] = to_numpy(std::move(c.ref_sig_start), {size_t(n + 1)});
        d["ref_sig"] = to_numpy(std::move(c.ref_sig), {S2});
        d["pair_sig_start"] = to_numpy(std::move(c.pair_sig_start), {size_t(q.n + 1)});
        d["pair_sig"] = to_numpy(std::move(c.pair_sig), {Q});
        d["pair_ref_sig_start"] = to_numpy(std::move(c.pair_ref_sig_start), {size_t(q.n + 1)});
        d["pair_ref_sig"] = to_numpy(std::move(c.pair_ref_sig), {Q2});
        return d;
      },
      nb::arg("xy").noconvert(), nb::arg("ring_start").noconvert(), nb::arg("shape_start").noconvert(),
      nb::arg("roles").noconvert().none(), nb::arg("W").noconvert(), nb::arg("H").noconvert(), nb::arg("G"),
      nb::arg("ref_bits").noconvert(), nb::arg("ref_hw").noconvert(), nb::arg("ref_start").noconvert(),
      nb::arg("ref_of").noconvert(), nb::arg("min_area").noconvert(), nb::arg("pairs").noconvert().none() = nb::none(),
      nb::arg("ref_pairs").noconvert().none() = nb::none(), nb::arg("pair_min_area").noconvert().none() = nb::none());
}
