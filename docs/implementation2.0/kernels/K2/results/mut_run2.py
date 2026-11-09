import subprocess, os, shutil, sys
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
muts = [
 ("const i64 v = 1 + sep(S(q), u);", "const i64 v = sep(S(q), u);"),
 ("const i64 src = low ? i64(32 * j + 31 - clz(low)) : last;", "const i64 src = low ? i64(32 * j + 31 - clz(low)) : -1;"),
 ("if (src >= 0 && src - x < o[x]) o[x] = i32(src - x);", "if (src >= 0 && src - x < o[x] && j > 0) o[x] = i32(src - x);"),
 ("while (q >= 0 && F(T(q), S(q)) > F(T(q), u)) --q;", "while (q >= 0 && F(T(q), S(q)) >= F(T(q), u)) --q;"),
]
src = open('k2score.inc').read()
for i, (a, b) in enumerate(muts, 1):
    assert a in src, a
    d = f'build/mut/n{i}'; os.makedirs(d + '/vendor', exist_ok=True)
    for f in os.listdir('vendor'): shutil.copy('vendor/' + f, d + '/vendor/')
    for f in ['k2.cpp', 'k2.h', 'k2ref.cpp', 'k2ref.h', 'k2test.cpp', 'k2gen.h']: shutil.copy(f, d)
    open(d + '/k2score.inc', 'w').write(src.replace(a, b, 1))
    r = subprocess.run(['clang++', '-std=c++17', '-O1', '-ffp-contract=off', d + '/k2.cpp', d + '/k2ref.cpp', d + '/k2test.cpp', '-o', d + '/t'], capture_output=True, text=True)
    if r.returncode: print(i, 'compile error', r.stderr[:300]); continue
    out = subprocess.run([d + '/t'], capture_output=True, text=True).stdout
    print(i, a[:50].replace('\n', ' '), '->', [l for l in out.splitlines() if 'failures' in l])
