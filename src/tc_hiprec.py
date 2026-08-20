"""Extended-precision (80-bit) two-chord E[A_4^2 & convex] for regular m-gons, two settings
each so the pair difference is an honest error bar."""
import json, os, sys, time
import numpy as np
import n6_twochord_polygon as TC

ms = [int(x) for x in sys.argv[1:]]
W = int(os.environ.get("W", "14"))
out = {}
for m in ms:
    V, q = TC.regular_mgon(m), TC.qsym_of(m)
    t0 = time.time()
    _, a = TC.two_chord_polygon(V, 40, 40, 6, q, workers=W, grade=1, dtname="ld")
    _, b = TC.two_chord_polygon(V, 52, 52, 6, q, workers=W, grade=1, dtname="ld")
    _, c2 = TC.two_chord_polygon(V, 44, 44, 6, q, workers=W, grade=2, dtname="ld")
    S = lambda z: repr(z).split("'")[1]
    print(f"m={m:3d}  A={S(a)}\n       B={S(b)}\n       C={S(c2)}  ({time.time()-t0:.0f}s)",
          flush=True)
    out[m] = dict(q=q, A=S(a), B=S(b), C=S(c2))
json.dump(out, open(os.path.join(os.path.dirname(__file__), "..", "results",
                                 "n6_tc_hiprec.json"), "w"), indent=1)
