"""Assemble and self-check the exact n=6 ingredient table for triangle, square, disk.

P_6 = 1 - 6 E[A_5] + 15 E[A_4^2] - 20 E[A_3^3],   E[A_4^2] = E[A_4^2 & convex] + 4 E[A_3^3].

Provenance of each ingredient:
  E[A_5]   : B-P width-function route (n6_bp_polygon.EA_hull_polygon); rational for triangle
             (43/180) and square (79/360) -- Buchta 1984 values; disk from E[N_6]=30 E[c^4].
  E[A_3^3] : B-P width-function route (n6_bp_polygon.EA3k_polygon), PSLQ; MC-confirmed.
             triangle 31/9000, square 137/72000, disk 1001/(6400 pi^4).
  E[A_4^2] : for triangle/square, the value forced by the decomposition + Valtr's exact P_6,
             independently CONFIRMED by Monte Carlo (mc_a4sq.py, rel.diff ~1e-4); disk from
             Marckert 2017.
Loop check: P_6 assembled == Valtr (triangle 91/900, square 49/400) / Marckert (disk).
"""
import sympy as sp
from fractions import Fraction

pi = sp.pi

DATA = {
    "triangle": dict(
        EA5=sp.Rational(43, 180), EA33=sp.Rational(31, 9000), EA42=sp.Rational(181, 4500),
        P6=sp.Rational(91, 900), P6_src="Valtr 1996"),
    "square": dict(
        EA5=sp.Rational(79, 360), EA33=sp.Rational(137, 72000), EA42=sp.Rational(859, 27000),
        P6=sp.Rational(49, 400), P6_src="Valtr 1995"),
    "disk": dict(
        EA5=7 * (2400 * pi**2 - 3289) / (6912 * pi**4),
        EA33=sp.Rational(1001) / (6400 * pi**4),
        EA42=(2400 * pi**2 + 31031) / (19200 * pi**4),
        P6=1 - (146400 * pi**2 - 473473) / (11520 * pi**4), P6_src="Marckert 2017"),
}


def report():
    rows = []
    allok = True
    for name, d in DATA.items():
        EA5, EA33, EA42, P6 = d["EA5"], d["EA33"], d["EA42"], d["P6"]
        P6_asm = sp.simplify(1 - 6 * EA5 + 15 * EA42 - 20 * EA33)
        ok = sp.simplify(P6_asm - P6) == 0
        allok &= ok
        EA42conv = sp.simplify(EA42 - 4 * EA33)
        print(f"=== {name} ===")
        print(f"  E[A_5]            = {EA5}   = {sp.N(EA5,12)}")
        print(f"  E[A_3^3]          = {EA33}   = {sp.N(EA33,12)}")
        print(f"  E[A_4^2]          = {EA42}   = {sp.N(EA42,12)}")
        print(f"  E[A_4^2 & convex] = {EA42conv}   = {sp.N(EA42conv,12)}")
        print(f"  P_6 assembled     = {sp.nsimplify(P6_asm)}   = {sp.N(P6,12)}")
        print(f"  P_6 ({d['P6_src']}) -> loop {'OK' if ok else 'FAIL'}")
        rows.append(dict(body=name, EA5=str(EA5), EA33=str(EA33), EA42=str(EA42),
                         EA42_convex=str(EA42conv), P6=str(P6), P6_decimal=float(sp.N(P6, 15)),
                         loop_ok=bool(ok), P6_source=d["P6_src"]))
    print("\nALL LOOPS OK" if allok else "\nLOOP FAILURE")
    return rows, allok


if __name__ == "__main__":
    import json, os
    rows, ok = report()
    out = os.path.join(os.path.dirname(__file__), "..", "results", "n6_ingredients.json")
    json.dump({"rows": rows, "all_loops_ok": ok,
               "decomposition": "P_6 = 1 - 6 E[A_5] + 15 E[A_4^2] - 20 E[A_3^3]",
               "note": "E[A_5], E[A_3^3] exact via Blaschke-Petkantschin width route; "
                       "E[A_4^2] Valtr/Marckert-implied, MC-confirmed (rel.diff ~1e-4)."},
              open(out, "w"), indent=1)
    print("wrote", os.path.normpath(out))
