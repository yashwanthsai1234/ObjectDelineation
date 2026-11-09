// K3 CPU reference (see k3_ref.h). Build with -ffp-contract=off (Makefile) so no a*b+c is fused implicitly.
#include "k3_ref.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <thread>

using namespace k3;

namespace k3ref {

PatchGeom geom(const PreParams& p) {
  PatchGeom g;
  g.H = p.H;
  g.W = p.W;
  g.ps = p.ps;
  g.ms = p.ms;
  g.tps = p.tps;
  g.gh = p.H / p.ps;
  g.gw = p.W / p.ps;
  g.F = p.C * p.tps * p.ps * p.ps;
  g.G = g.gh * g.gw;
  if (p.H % (p.ps * p.ms) || p.W % (p.ps * p.ms)) throw std::invalid_argument("H, W must be multiples of ps*ms");
  return g;
}

static inline uint8_t px_at(const uint8_t* img, const PreParams& p, int c, int y, int x) {
  return p.hwc ? img[((size_t)y * p.W + x) * p.C + c] : img[((size_t)c * p.H + y) * p.W + x];
}

void preprocess_ref(const uint8_t* src, const int* idx, int n_out, const PreParams& p, float* out) {
  const PatchGeom g = geom(p);
  const size_t img_bytes = (size_t)p.H * p.W * p.C;
  std::vector<float> chw((size_t)p.C * p.H * p.W);
  size_t o = 0;
  for (int n = 0; n < n_out; ++n) {
    const uint8_t* img = src + (size_t)(idx ? idx[n] : n) * img_bytes;
    // rescale_and_normalize on the whole image (HF: images.to(float32), sub_(mean'), div_(std'))
    for (int c = 0; c < p.C; ++c)
      for (int y = 0; y < p.H; ++y)
        for (int x = 0; x < p.W; ++x) {
          float v = (float)px_at(img, p, c, y, x);
          v = v - p.mean[c];
          v = v / p.std_[c];
          chw[((size_t)c * p.H + y) * p.W + x] = v;
        }
    // patchify: reshape (C, GH, MH, PY, GW, MW, PX) -> permute (GH, GW, MH, MW, C, PY, PX) -> expand T after C
    const int GH = g.gh / p.ms, GW = g.gw / p.ms;
    for (int gh = 0; gh < GH; ++gh)
      for (int gw = 0; gw < GW; ++gw)
        for (int mh = 0; mh < p.ms; ++mh)
          for (int mw = 0; mw < p.ms; ++mw)
            for (int c = 0; c < p.C; ++c)
              for (int t = 0; t < p.tps; ++t)
                for (int py = 0; py < p.ps; ++py)
                  for (int px = 0; px < p.ps; ++px) {
                    int y = (gh * p.ms + mh) * p.ps + py;
                    int x = (gw * p.ms + mw) * p.ps + px;
                    out[o++] = chw[((size_t)c * p.H + y) * p.W + x];
                  }
  }
}

void preprocess_emul(const uint8_t* src, const int* idx, int n_out, const PreParams& p, float* out) {
  const PatchGeom g = geom(p);
  const size_t img_bytes = (size_t)p.H * p.W * p.C;
  for (int n = 0; n < n_out; ++n) {
    const uint8_t* img = src + (size_t)(idx ? idx[n] : n) * img_bytes;
    for (int q = 0; q < g.G; ++q)
      for (int f = 0; f < g.F; ++f) {
        int c, y, x;
        k3_patch_src(g, q, f, &c, &y, &x);
        out[((size_t)n * g.G + q) * g.F + f] = k3_norm(px_at(img, p, c, y, x), p.mean[c], p.std_[c]);
      }
  }
}

// ---------------- rope ----------------
// One sequence of real tokens -> positions (3 x L) and the final current_pos. HF's loop, literally.
static void rope_core(const int64_t* types, int L, const std::vector<Grid>& grids, int merge, int64_t* pos,
                      int64_t stride, int64_t* end) {
  size_t gi = 0;
  int64_t current = 0;
  int k = 0;
  while (k < L) {
    int64_t ty = types[k];
    int e = k;
    while (e < L && types[e] == ty) ++e;
    if (ty == 0) {
      for (int j = k; j < e; ++j)
        for (int a = 0; a < 3; ++a) pos[a * stride + j] = (j - k) + current;
      current += e - k;
    } else if (ty == 1) {
      if (gi >= grids.size()) throw std::invalid_argument("more image groups than grids");
      const Grid g = grids[gi++];
      const int64_t lt = g[0], lh = g[1] / merge, lw = g[2] / merge;
      if (lt * lh * lw != e - k) throw std::invalid_argument("image group length != grid tokens");
      int j = k;
      for (int64_t t = 0; t < lt; ++t)
        for (int64_t h = 0; h < lh; ++h)
          for (int64_t w = 0; w < lw; ++w, ++j) {
            pos[0 * stride + j] = t + current;  // (arange(t) * 1.0).long() + start
            pos[1 * stride + j] = h + current;
            pos[2 * stride + j] = w + current;
          }
      current += std::max(g[1], g[2]) / merge;
    } else {
      throw std::invalid_argument("video tokens are not supported");
    }
    k = e;
  }
  *end = current;
}

int64_t rope_ref_row(const int64_t* ids, const int64_t* mask, const int64_t* types, int W,
                     const std::vector<Grid>& grids, int merge, int64_t* pos) {
  (void)ids;
  std::vector<int64_t> rtypes;
  std::vector<int> cols;
  for (int j = 0; j < W; ++j)
    if (mask[j]) {
      rtypes.push_back(types[j]);
      cols.push_back(j);
    }
  const int L = (int)cols.size();
  std::vector<int64_t> rp(3 * (size_t)std::max(L, 1));
  int64_t end = 0;
  rope_core(rtypes.data(), L, grids, merge, rp.data(), L, &end);
  for (int a = 0; a < 3; ++a)
    for (int j = 0; j < W; ++j) pos[a * W + j] = 0;
  int64_t mx = 0;
  for (int a = 0; a < 3; ++a)
    for (int i = 0; i < L; ++i) {
      pos[a * W + cols[i]] = rp[a * L + i];
      mx = std::max(mx, rp[a * L + i]);
    }
  return mx + 1 - L;
}

RopeTemplate rope_template(const int64_t* ids, const int64_t* types, int Lp, const std::vector<Grid>& grids,
                           int merge) {
  RopeTemplate t;
  t.Lp = Lp;
  t.ids.assign(ids, ids + Lp);
  t.types.assign(types, types + Lp);
  t.pos.assign(3 * (size_t)Lp, 0);
  rope_core(types, Lp, grids, merge, t.pos.data(), Lp, &t.end);
  return t;
}

void rope_emul(int B, int W, const int64_t* ids, const int64_t* mask, const int64_t* types, const RopeTemplate& t,
               int64_t* pos, int* err) {
  *err = 0;
  for (int b = 0; b < B; ++b) {
    int64_t r = 0;
    for (int j = 0; j < W; ++j) {
      const size_t o = (size_t)b * W + j;
      if (!mask[o]) {
        for (int a = 0; a < 3; ++a) pos[(size_t)a * B * W + o] = 0;
        continue;
      }
      if (r < t.Lp) {
        if (ids[o] != t.ids[r] || types[o] != t.types[r]) *err |= 1;
      } else if (types[o] != 0) {
        *err |= 2;
      }
      for (int a = 0; a < 3; ++a) pos[(size_t)a * B * W + o] = k3_rope_pos(r, t.Lp, t.pos.data() + a * t.Lp, t.end);
      ++r;
    }
    if (r < t.Lp) *err |= 4;
  }
}

// ---------------- fused log-probabilities ----------------
int lse_tiles(int V) { return (V + LSE_BN - 1) / LSE_BN; }
int lse_tiles_per_split(int V, int S) { return (lse_tiles(V) + S - 1) / S; }

static void lse_rows(const uint16_t* H, const uint16_t* Wv, const int32_t* tgt, int n0, int n1, int V, int D,
                     float inv_t, int S, float* logp, float* ent) {
  const int NT = lse_tiles(V), TPS = lse_tiles_per_split(V, S);
  std::vector<float> z(V), hf(D), ms(S), ss(S), xs(S);
  for (int n = n0; n < n1; ++n) {
    for (int k = 0; k < D; ++k) hf[k] = bf16_to_f32(H[(size_t)n * D + k]);
    for (int v = 0; v < V; ++v) {
      const uint16_t* w = Wv + (size_t)v * D;
      float acc = 0.0f;
      for (int k = 0; k < D; ++k) acc = K3_FMA(hf[k], bf16_to_f32(w[k]), acc);
      z[v] = K3_MUL(acc, inv_t);
    }
    for (int s = 0; s < S; ++s) {
      LseState st = lse_init();
      for (int t = s * TPS; t < std::min(NT, (s + 1) * TPS); ++t) {
        const int v0 = t * LSE_BN, nv = std::min(LSE_BN, V - v0);
        float mq[LSE_LANES], E[LSE_LANES], X[LSE_LANES];
        for (int q = 0; q < LSE_LANES; ++q) {
          mq[q] = -K3_INF;
          for (int i = 0; i < LSE_PER_LANE; ++i) {
            int j = LSE_LANES * i + q;
            if (j < nv) mq[q] = k3_max(mq[q], z[v0 + j]);
          }
        }
        const float Mt = k3_max(k3_max(mq[0], mq[1]), k3_max(mq[2], mq[3]));
        const float new_m = k3_max(st.m, Mt);
        for (int q = 0; q < LSE_LANES; ++q) {
          E[q] = 0.0f;
          X[q] = 0.0f;
          if (new_m == -K3_INF) continue;
          for (int i = 0; i < LSE_PER_LANE; ++i) {
            int j = LSE_LANES * i + q;
            if (j < nv) lse_lane_accumulate(z[v0 + j], new_m, E[q], X[q]);
          }
        }
        lse_fold(st, new_m, K3_ADD(K3_ADD(E[0], E[1]), K3_ADD(E[2], E[3])),
                 K3_ADD(K3_ADD(X[0], X[1]), K3_ADD(X[2], X[3])));
      }
      ms[s] = st.m;
      ss[s] = st.s;
      xs[s] = st.x;
    }
    const bool valid = tgt[n] >= 0 && tgt[n] < V;
    LseOut o = lse_finish(ms.data(), ss.data(), xs.data(), S, 1, valid ? z[tgt[n]] : 0.0f, valid);
    logp[n] = o.logp;
    if (ent) ent[n] = o.entropy;
  }
}

template <class F>
static void parallel_rows(int N, F f) {
  int nt = (int)std::max(1u, std::thread::hardware_concurrency());
  nt = std::min(nt, std::max(N, 1));
  std::vector<std::thread> th;
  const int per = (N + nt - 1) / nt;
  for (int i = 0; i < nt; ++i) {
    int a = i * per, b = std::min(N, a + per);
    if (a < b) th.emplace_back([=] { f(a, b); });
  }
  for (auto& t : th) t.join();
}

void lse_ref(const uint16_t* H, const uint16_t* Wv, const int32_t* tgt, int N, int V, int D, float inv_t, int S,
             float* logp, float* ent) {
  parallel_rows(N, [&](int a, int b) { lse_rows(H, Wv, tgt, a, b, V, D, inv_t, S, logp, ent); });
}

void lse_ref_f64(const uint16_t* H, const uint16_t* Wv, const int32_t* tgt, int N, int V, int D, double inv_t,
                 double* logp, double* ent) {
  parallel_rows(N, [&](int a, int b) {
    std::vector<double> z(V);
    for (int n = a; n < b; ++n) {
      double m = -INFINITY;
      for (int v = 0; v < V; ++v) {
        double acc = 0;
        for (int k = 0; k < D; ++k)
          acc += (double)bf16_to_f32(H[(size_t)n * D + k]) * (double)bf16_to_f32(Wv[(size_t)v * D + k]);
        z[v] = acc * inv_t;
        m = std::max(m, z[v]);
      }
      double s = 0, x = 0;
      for (int v = 0; v < V; ++v) {
        double e = std::exp(z[v] - m);
        s += e;
        x += e * z[v];
      }
      const bool valid = tgt[n] >= 0 && tgt[n] < V;
      logp[n] = valid ? z[tgt[n]] - m - std::log(s) : 0.0;
      if (ent) ent[n] = m + std::log(s) - x / s;
    }
  });
}

// ---------------- IS ratio ----------------
static float tree(std::vector<float>& a, int op) {
  for (int stride = IS_THREADS / 2; stride > 0; stride >>= 1)
    for (int i = 0; i < stride; ++i) a[i] = red_apply(op, a[i], a[i + stride]);
  return a[0];
}

static const int STAT_OP[IS_NSTATS] = {RED_SUM, RED_SUM, RED_MAX, RED_SUM, RED_SUM, RED_MIN, RED_MAX};
static const float STAT_INIT[IS_NSTATS] = {0.0f, 0.0f, -K3_INF, 0.0f, 0.0f, K3_INF, -K3_INF};

void is_ref(const float* lp, const float* sp, const uint8_t* mask, int B, int T, int mode, float lo, float hi,
            float* ratio, float* row_stats, float* stats) {
  const bool seq = mode == SEQUENCE_TRUNCATE || mode == SEQUENCE_MASK;
  std::vector<float> part[IS_NSTATS + 1];
  for (auto& v : part) v.resize(IS_THREADS);
  for (int b = 0; b < B; ++b) {
    for (int tid = 0; tid < IS_THREADS; ++tid) {
      float a[IS_NSTATS];
      for (int k = 0; k < IS_NSTATS; ++k) a[k] = STAT_INIT[k];
      float dacc = 0.0f;
      for (int j = tid; j < T; j += IS_THREADS) {
        const size_t o = (size_t)b * T + j;
        const float d = is_diff(lp[o], sp[o], mask[o]);
        if (mask[o] && !k3_isnan(sp[o])) {
          float ad = K3_SUB(lp[o], sp[o]);
          ad = ad < 0.0f ? -ad : ad;
          a[0] = K3_ADD(a[0], ad);
          a[1] = K3_ADD(a[1], 1.0f);
          a[2] = k3_max(a[2], ad);
        }
        if (seq) {
          dacc = K3_ADD(dacc, d);
        } else {
          const float r = is_bound(k3_expf(d), mode, lo, hi);
          ratio[o] = r;
          if (mask[o]) {
            a[3] = K3_ADD(a[3], r);
            a[4] = K3_ADD(a[4], 1.0f);
            a[5] = k3_min(a[5], r);
            a[6] = k3_max(a[6], r);
          }
        }
      }
      for (int k = 0; k < IS_NSTATS; ++k) part[k][tid] = a[k];
      part[IS_NSTATS][tid] = dacc;
    }
    float* rs = row_stats + (size_t)b * IS_NSTATS;
    for (int k = 0; k < 3; ++k) rs[k] = tree(part[k], STAT_OP[k]);
    if (seq) {
      const float D = tree(part[IS_NSTATS], RED_SUM);
      const float r = is_bound(k3_expf(D), mode, lo, hi);
      ratio[b] = r;
      rs[3] = r;
      rs[4] = 1.0f;
      rs[5] = r;
      rs[6] = r;
    } else {
      for (int k = 3; k < IS_NSTATS; ++k) rs[k] = tree(part[k], STAT_OP[k]);
    }
  }
  // over rows: thread tid folds rows tid, tid + 256, ... in order, then the tree
  for (int k = 0; k < IS_NSTATS; ++k) {
    for (int tid = 0; tid < IS_THREADS; ++tid) {
      float a = STAT_INIT[k];
      for (int b = tid; b < B; b += IS_THREADS) a = red_apply(STAT_OP[k], a, row_stats[(size_t)b * IS_NSTATS + k]);
      part[k][tid] = a;
    }
    stats[k] = tree(part[k], STAT_OP[k]);
  }
}

}  // namespace k3ref
