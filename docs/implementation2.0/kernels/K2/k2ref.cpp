// k2ref.cpp -- see k2ref.h. Plain loops, no bit tricks, no shared code with k2score.inc.
#include "k2ref.h"

#include <algorithm>
#include <stdexcept>

namespace k2ref {

std::vector<Mask> unpack(const k2::Bits& b) {
  std::vector<Mask> out(size_t(b.n));
  for (i64 m = 0; m < b.n; ++m) {
    Mask& k = out[size_t(m)];
    k.H = b.H[m], k.W = b.W[m];
    k.px.assign(size_t(k.H) * k.W, 0);
    const i32 wpr = (k.W + 31) / 32;
    for (i32 y = 0; y < k.H; ++y)
      for (i32 x = 0; x < k.W; ++x)
        k.px[size_t(y) * k.W + x] = (b.words[b.word_start[m] + i64(y) * wpr + x / 32] >> (x % 32)) & 1u;
  }
  return out;
}

void pack(const std::vector<Mask>& ms, std::vector<u32>& words, std::vector<i64>& word_start, std::vector<i32>& H,
          std::vector<i32>& W) {
  words.clear(), word_start.assign(1, 0), H.clear(), W.clear();
  for (const Mask& m : ms) {
    const i32 wpr = (m.W + 31) / 32;
    const size_t at = words.size();
    words.resize(at + size_t(m.H) * wpr, 0u);
    for (i32 y = 0; y < m.H; ++y)
      for (i32 x = 0; x < m.W; ++x)
        if (m.px[size_t(y) * m.W + x]) words[at + size_t(y) * wpr + x / 32] |= 1u << (x % 32);
    word_start.push_back(i64(words.size()));
    H.push_back(m.H), W.push_back(m.W);
  }
}

Mask boundary(const Mask& m) {
  Mask b = m;
  for (i32 y = 0; y < m.H; ++y)
    for (i32 x = 0; x < m.W; ++x) {
      const bool f = m.at(x, y);
      const bool inner = m.at(x - 1, y) && m.at(x + 1, y) && m.at(x, y - 1) && m.at(x, y + 1);
      b.px[size_t(y) * m.W + x] = f && !inner;
    }
  return b;
}

std::vector<i32> edt_brute(const Mask& src) {
  std::vector<std::pair<i32, i32>> pts;
  for (i32 y = 0; y < src.H; ++y)
    for (i32 x = 0; x < src.W; ++x)
      if (src.at(x, y)) pts.push_back({x, y});
  std::vector<i32> d(size_t(src.H) * src.W, k2::kNoSource);
  if (pts.empty()) return d;
  for (i32 y = 0; y < src.H; ++y)
    for (i32 x = 0; x < src.W; ++x) {
      i64 best = -1;
      for (const auto& p : pts) {
        const i64 dx = x - p.first, dy = y - p.second, v = dx * dx + dy * dy;
        if (best < 0 || v < best) best = v;
      }
      d[size_t(y) * src.W + x] = i32(best);
    }
  return d;
}

Directed directed(const Mask& from, const Mask& to, i64 tau2, i64 p_num, i64 p_den) {
  const Mask bf = boundary(from), bt = boundary(to);
  const std::vector<i32> dist = edt_brute(bt);
  std::vector<i64> v;
  Directed r;
  for (i32 y = 0; y < bf.H; ++y)
    for (i32 x = 0; x < bf.W; ++x)
      if (bf.at(x, y)) {
        const i64 d = dist[size_t(y) * bf.W + x];
        v.push_back(d);
        r.hit += d <= tau2;
      }
  r.n = i64(v.size());
  if (v.empty()) return r;
  std::sort(v.begin(), v.end());
  r.rank = p_num * (r.n - 1) / p_den;
  r.lo = v[size_t(r.rank)];
  r.hi = v[size_t(std::min(r.rank + 1, r.n - 1))];
  return r;
}

static double iou_of(i64 inter, i64 a, i64 b) {
  const i64 u = a + b - inter;
  return u ? double(inter) / double(u) : 1.0;
}

Set set_stats(const std::vector<Mask>& masks, i32 K) {
  const i64 M = i64(masks.size()), k = M - K;
  if (K < 1 || k < 1) throw std::invalid_argument("K >= 1 and k >= 1");
  Set s;
  s.inter.assign(size_t(M * M), 0);
  for (i64 i = 0; i < M; ++i)
    for (i64 j = 0; j < M; ++j) {
      i64 c = 0;
      for (size_t p = 0; p < masks[size_t(i)].px.size(); ++p) c += masks[size_t(i)].px[p] && masks[size_t(j)].px[p];
      s.inter[size_t(i * M + j)] = c;
    }
  auto I = [&](i64 i, i64 j) {
    return iou_of(s.inter[size_t(i * M + j)], s.inter[size_t(i * M + i)], s.inter[size_t(j * M + j)]);
  };
  double sxy = 0, sxx = 0, syy = 0;
  for (i64 i = 0; i < K; ++i)
    for (i64 j = K; j < M; ++j) sxy += 1.0 - I(i, j);
  for (i64 i = 0; i < K; ++i)
    for (i64 j = 0; j < K; ++j) sxx += 1.0 - I(i, j);
  for (i64 i = K; i < M; ++i)
    for (i64 j = K; j < M; ++j) syy += 1.0 - I(i, j);
  const double g = 2.0 * (sxy / double(K * k)) - sxx / double(K * K) - syy / double(k * k);
  s.ged = g > 0.0 ? g : 0.0;
  // all k^K assignments in lexicographic order; the first strict maximum of the fixed-point total
  const i64 base = K / k, extra = K % k;
  auto Q = [&](i64 i, i64 r) -> u64 {
    const i64 in = s.inter[size_t(i * M + K + r)], a = s.inter[size_t(i * M + i)],
              b = s.inter[size_t((K + r) * M + K + r)], u = a + b - in;
    return u ? (u64(in) << 32) / u64(u) : (u64(1) << 32);
  };
  std::vector<i32> cur(size_t(K), 0), cnt(size_t(k), 0);
  bool found = false;
  for (;;) {
    std::fill(cnt.begin(), cnt.end(), 0);
    for (i32 c : cur) ++cnt[size_t(c)];
    i64 plus = 0;
    bool ok = true;
    for (i32 c : cnt) {
      if (c < base || c > base + (extra ? 1 : 0)) ok = false;
      plus += c == base + 1;
    }
    if (ok && (extra == 0 || plus == extra)) {
      u64 tot = 0;
      for (i64 i = 0; i < K; ++i) tot += Q(i, cur[size_t(i)]);
      if (!found || tot > s.hm_q) s.hm_q = tot, s.choice = cur, found = true;
    }
    i64 p = K - 1;  // next sequence
    while (p >= 0 && cur[size_t(p)] == k - 1) cur[size_t(p--)] = 0;
    if (p < 0) break;
    ++cur[size_t(p)];
  }
  double sum = 0;
  for (i64 i = 0; i < K; ++i) sum += I(i, K + s.choice[size_t(i)]);
  s.hm_iou = sum / double(K);
  return s;
}

}  // namespace k2ref
