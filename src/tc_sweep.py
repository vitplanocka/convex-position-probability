"""Two-chord E[A_4^2 & convex] for regular m-gons over a wide m range, two grids each."""
import json, os, sys, time
import n6_twochord_polygon as TC
ms = [int(x) for x in sys.argv[1:]]
out = {}
for m in ms:
    V, q = TC.regular_mgon(m), TC.qsym_of(m)
    t0 = time.time()
    _, a = TC.two_chord_polygon(V, 28, 28, 6, q, workers=5)
    _, b = TC.two_chord_polygon(V, 36, 36, 6, q, workers=5)
    print(f"m={m:3d} q={q:2d}  {float(b):.17g}   spread={abs(float(a)-float(b)):.2e}  "
          f"({time.time()-t0:.0f}s)", flush=True)
    out[m] = dict(q=q, EA42_convex=float(b), spread=abs(float(a) - float(b)))
json.dump(out, open(os.path.join(os.path.dirname(__file__), "..", "results",
                                 "n6_tc_sweep.json"), "w"), indent=1)
