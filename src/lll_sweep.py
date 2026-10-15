"""Parameter sweep for lll_search.py; logs every run and keeps the best witness per (field, sign).
usage: python src/lll_sweep.py OUTDIR"""
import itertools, json, re, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
ZF = {'Q(sqrt-3)': 'data/zeros/hp_Qsqrt-3_T1000_D100.json', 'Q(sqrt5)': 'data/zeros/hp_Qsqrt5_T1000_D100.json'}
grid = list(itertools.product(['Q(sqrt5)', 'Q(sqrt-3)'], [1, -1], [600, 1000], [20, 30, 40], [5, 7], ['c', 'sqrtc']))
best = {}
log = open(out / 'sweep.log', 'a', encoding='utf-8')
t0 = time.time()
for i, (field, sg, T, n, f, wt) in enumerate(grid):
    nu = n * f
    cmd = [sys.executable, 'src/lll_search.py', field, str(T), str(n), str(nu), '--sign', str(sg), '--zeros', ZF[field],
           '--weight', wt]
    r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    m = re.search(r'BEST: .*h\*\((-?\d+)\) = (\S+)', r.stdout)
    line = f"{field} sign={sg:+d} T={T} n={n} nu={nu} w={wt}: " + (f"h*={m.group(2)} t={m.group(1)}" if m else "no candidate")
    log.write(line + '\n'); log.flush()
    if m:
        v = float(m.group(2))
        key = f"{field}|{sg:+d}"
        if key not in best or sg * v > sg * best[key]['h']:
            best[key] = {'h': v, 't': m.group(1), 'T': T, 'n': n, 'nu': nu, 'weight': wt}
            json.dump(best, open(out / 'best.json', 'w'), indent=1)
    print(f"[{i + 1}/{len(grid)} {time.time() - t0:.0f}s] {line[:110]}", flush=True)
print(json.dumps(best, indent=1))
