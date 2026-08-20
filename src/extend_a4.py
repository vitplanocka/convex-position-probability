"""Extend the n=6 m-gon table to larger m to sharpen the large-m asymptotics
P_6(disk) - P_6(m) = a4/m^4 + a6/m^6 + ...   (for P_5 the exact answer is a4 = 7 pi^2/18.)"""
import json, os, sys, time
import mpmath as mp
import n6_bp_polygon as BP
import n6_twochord_polygon as TC
mp.mp.dps = 40

def one(m):
    t0 = time.time()
    V, area = BP.regular_mgon(m)
    EA33 = BP.EA3k_polygon(V, 3, area)
    EA5 = BP.EA_hull_polygon(V, area, 5)[5]
    Vf, q = TC.regular_mgon(m), TC.qsym_of(m)
    _, a = TC.two_chord_polygon(Vf, 36, 36, 6, q, workers=int(os.environ.get("TCW","3")), grade=1, dtname="ld")
    _, b = TC.two_chord_polygon(Vf, 44, 44, 6, q, workers=int(os.environ.get("TCW","3")), grade=1, dtname="ld")
    S = lambda z: repr(z).split("'")[1]
    ea42 = mp.mpf(S(b))
    P6 = 1 - 6 * EA5 + 15 * ea42 + 40 * EA33
    return dict(m=m, EA5=mp.nstr(EA5, 25), EA33=mp.nstr(EA33, 25),
                EA42_convex=S(b), EA42_convex_alt=S(a), P6=mp.nstr(P6, 22),
                secs=round(time.time() - t0, 1))

if __name__ == "__main__":
    # one m per top-level process (nested pools are not allowed); merge the JSONs afterwards
    for m in (int(x) for x in sys.argv[1:]):
        r = one(m)
        print(json.dumps(r), flush=True)
        json.dump(r, open(f"../results/n6_mgon_large_m_{m}.json", "w"), indent=1)
