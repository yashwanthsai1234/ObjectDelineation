#!/bin/bash
# mutation check: each mutant of k2score.inc must make k2test fail
cd "$(dirname "$0")/../.."
i=0
while IFS='|' read -r from to; do
  i=$((i+1)); d=build/mut/m$i; mkdir -p $d/vendor; cp vendor/* $d/vendor/; cp k2.cpp k2.h k2ref.* k2test.cpp k2gen.h $d/
  python3 -c "
import sys; s=open('k2score.inc').read(); a,b=sys.argv[1],sys.argv[2]
assert a in s, a; open('$d/k2score.inc','w').write(s.replace(a,b,1))" "$from" "$to" || { echo "m$i: pattern missing"; continue; }
  clang++ -std=c++17 -O1 -ffp-contract=off $d/k2.cpp $d/k2ref.cpp $d/k2test.cpp -o $d/t 2>/dev/null || { echo "m$i: compile error"; continue; }
  out=$($d/t | grep -E "failures")
  echo "m$i [$from]: $out"
done <<'M'
!= ccx || grid_cell|!= ccx && grid_cell
const bool first = v < u;|const bool first = v > u;
if (wany(on)) return 0;|if (false && wany(on)) return 0;
if (wany(hit)) return true;|if (base == 0 && wany(hit)) return true;
at(h, kHoleOutsideOuter) = 1;|(void)0;
up & dn & left & right|dn & left & right
(j + 1 < wpr ? row[j + 1] << 31 : 0u)|0u
const i64 k = p_num * (n - 1) / p_den;|const i64 k = p_num * n / p_den;
(k + 1 >= n || le >= k + 2) ? ans : above|(k + 1 >= n) ? ans : above
if (c >= cap || B[st + rad] == kNeg || Q(i, r) + B[st + rad] != B[st]) continue;|if (c >= cap || B[st + rad] == kNeg || Q(i, r) + B[st + rad] < B[st]) continue;
if (b == kNeg || val > b) b = val;|if (b == kNeg || val < b) b = val;
B[st] = (sum == K && full) ? 0 : kNeg;|B[st] = (sum == K) ? 0 : kNeg;
if (below + 1 < o[i64(y) * W]) o[i64(y) * W] = below + 1;|(void)below;
const i64 w = 1 + sep(sr[q], u);|const i64 w = sep(sr[q], u);
h += dv <= tau2;|h += dv < tau2;
if (q / R < q % R) v += m[|if (q / R <= q % R) v += m[
bad += popc(wballot(nz));|bad += nz && a > 0;
M
