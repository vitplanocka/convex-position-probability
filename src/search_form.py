"""Search for a closed form X(m) = sum over groups G(a,b,J) = {c^j / (m^a s^b): j<=J}
for the regular-m-gon n=6 ingredients.  Solve on the first n points, check the rest."""
import itertools, json, sys
import mpmath as mp
mp.mp.dps = 40

J = json.load(open("../results/n6_twochord_polygon.json"))
ROWS = {r["m"]: r for r in J["rows"]}
try:
    SW = json.load(open("../results/n6_tc_sweep.json"))
except Exception:
    SW = {}

def cs(m):
    w = 2*mp.pi/m
    return mp.cos(w), mp.sin(w)

def make_basis(groups):
    fs = []
    for (a, b, d) in groups:
        for j in range(d+1):
            fs.append((a, b, j))
    return fs

def ev(fs, m):
    c, s = cs(m)
    return [c**j / (mp.mpf(m)**a * s**b) for (a, b, j) in fs]

def test(ms, ys, groups, tol):
    fs = make_basis(groups)
    n = len(fs)
    if n < 1 or n > len(ms) - 2:
        return None
    A = mp.matrix([ev(fs, ms[i]) for i in range(n)])
    b = mp.matrix([ys[i] for i in range(n)])
    try:
        x = mp.lu_solve(A, b)
    except Exception:
        return None
    res = mp.mpf(0)
    for i in range(len(ms)):
        r = mp.fdot(ev(fs, ms[i]), x) - ys[i]
        res = max(res, abs(r)/abs(ys[i]))
    return (res, x, fs) if res < tol else None

AB = [(2,2),(4,4),(6,6),(3,3),(5,5),(4,2),(2,4),(6,4),(4,6),(3,2),(2,3),(5,4),(4,5),
      (6,5),(5,6),(3,4),(4,3),(6,3),(3,6),(8,8),(6,2),(2,6),(8,6),(6,8)]

def run(name, ms, ys, tol, maxpar):
    print(f"### {name}: {len(ms)} points, tol {mp.nstr(tol,2)}, <= {maxpar} params")
    found = []
    for ng in (1, 2, 3):
        for combo in itertools.combinations(AB, ng):
            for ds in itertools.product(range(0, 6), repeat=ng):
                groups = [(a, b, d) for (a, b), d in zip(combo, ds)]
                npar = sum(d+1 for _, _, d in groups)
                if npar > maxpar or npar > len(ms)-2:
                    continue
                r = test(ms, ys, groups, tol)
                if r:
                    found.append((npar, float(r[0]), groups, r[1]))
        if found:
            break
    found.sort()
    for npar, res, groups, x in found[:4]:
        print(f"   FOUND {groups}  npar={npar} maxrelres={res:.2e}")
        print("      coeffs:", [mp.nstr(v, 20) for v in x])
    if not found:
        print("   nothing")
    return found

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    maxpar = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    ms = sorted(ROWS)
    if which in ("all", "EA33"):
        run("E[A_3^3]", ms, [mp.mpf(ROWS[m]["EA33"]) for m in ms], mp.mpf(10)**-18, maxpar)
    if which in ("all", "EA5"):
        run("E[A_5]", ms, [mp.mpf(ROWS[m]["EA5"]) for m in ms], mp.mpf(10)**-18, maxpar)
    if which in ("all", "EA42"):
        mw = sorted(int(k) for k in SW) or ms
        yw = [mp.mpf(repr(SW[str(m)]["EA42_convex"])) for m in mw] if SW else \
             [mp.mpf(ROWS[m]["EA42_convex"]) for m in mw]
        run("E[A_4^2&conv]", mw, yw, mp.mpf(10)**-13, maxpar)
