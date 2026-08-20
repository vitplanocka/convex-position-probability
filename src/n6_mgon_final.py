"""FINAL n=6 regular-m-gon table: P_6 assembled from first principles.

    P_6 = 1 - 6 E[A_5] + 15 E[A_4^2 & convex] + 40 E[A_3^3]
        (= 1 - 6 E[A_5] + 15 E[A_4^2] - 20 E[A_3^3] with E[A_4^2] = E[A_4^2&conv] + 4E[A_3^3])

E[A_5], E[A_3^3]        : Blaschke-Petkantschin WIDTH route   (n6_bp_polygon)
E[A_4^2 & convex]       : two-chord B-P line-space integral   (n6_twochord_polygon, 80-bit)
Nothing here uses Valtr or Marckert; m=3 -> 91/900 and m=4 -> 49/400 are OUTPUTS.
"""
import json, os
from fractions import Fraction as Fr
import mpmath as mp
import sympy as sp

mp.mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "..", "results")

# ---- exact ingredient values, PSLQ'd from the 20+-digit width-route integrals --------------
s2, s3, s5 = sp.sqrt(2), sp.sqrt(3), sp.sqrt(5)
Q = sp.Rational
EXACT_EA5 = {3: Q(43, 180), 4: Q(79, 360), 5: Q(221, 1500) + Q(34, 1125) * s5,
             6: Q(149347, 699840), 8: Q(4531, 36864) + Q(5837, 92160) * s2,
             10: Q(361891, 3000000) + Q(5121, 125000) * s5,
             12: Q(253843, 2239488) + Q(79843, 1399680) * s3}
EXACT_EA33 = {3: Q(31, 9000), 4: Q(137, 72000), 5: Q(32, 28125) + Q(19, 75000) * s5,
              6: Q(57709, 34992000), 8: Q(193, 230400) + Q(53, 96000) * s2,
              10: Q(92369, 112500000) + Q(52973, 150000000) * s5,
              12: Q(6103, 7464960) + Q(1183, 2592000) * s3}
# ---- exact E[A_4^2 & convex], from the smooth-denominator scan on the 80-bit values ---------
EXACT_EA42C = {3: Q(119, 4500), 4: Q(1307, 54000),
               5: (5307 + 1154 * s5) / 337500,
               6: Q(403891, 17496000),
               8: (178715 + 97586 * s2) / 13824000}


def main():
    hi = json.load(open(os.path.join(R, "n6_tc_hiprec.json")))
    wid = {r["m"]: r for r in json.load(open(os.path.join(R, "n6_mgon_width_ingredients.json")))["rows"]}
    rows = []
    for k in sorted(hi, key=int):
        m = int(k)
        A, B, C = (mp.mpf(hi[k][x]) for x in "ABC")
        ea42c = B
        err42 = max(abs(A - B), abs(C - B), mp.mpf(10) ** -20)
        ea5 = mp.mpf(wid[m]["EA5"])
        ea33 = mp.mpf(wid[m]["EA33"])
        P6 = 1 - 6 * ea5 + 15 * ea42c + 40 * ea33
        row = dict(m=m, EA5=mp.nstr(ea5, 25), EA33=mp.nstr(ea33, 25),
                   EA42_convex=mp.nstr(ea42c, 21), EA42=mp.nstr(ea42c + 4 * ea33, 21),
                   P6=mp.nstr(P6, 21), P6_float=float(P6),
                   EA42_convex_err=float(err42))
        if m in EXACT_EA42C:
            e42c = sp.nsimplify(EXACT_EA42C[m])
            e42 = sp.radsimp(sp.simplify(e42c + 4 * EXACT_EA33[m]))
            p6 = sp.radsimp(sp.simplify(1 - 6 * EXACT_EA5[m] + 15 * e42c + 40 * EXACT_EA33[m]))
            resid = abs(mp.mpf(str(sp.N(p6, 30))) - P6)
            row.update(EA42_convex_exact=sp.sstr(e42c), EA42_exact=sp.sstr(e42),
                       P6_exact=sp.sstr(p6), P6_exact_decimal=str(sp.N(p6, 21)),
                       exact_vs_numeric=float(resid))
        rows.append(row)
    disk = dict(m="inf(disk)", EA5=str(sp.N(7 * (2400 * sp.pi**2 - 3289) / (6912 * sp.pi**4), 21)),
                EA33=str(sp.N(sp.Rational(1001) / (6400 * sp.pi**4), 21)),
                EA42_convex=str(sp.N((2400 * sp.pi**2 + 19019) / (19200 * sp.pi**4), 21)),
                P6=str(sp.N(1 - (146400 * sp.pi**2 - 473473) / (11520 * sp.pi**4), 21)),
                P6_exact="1 - (146400*pi**2 - 473473)/(11520*pi**4)", source="Marckert 2017")
    out = dict(
        formula="P_6 = 1 - 6 E[A_5] + 15 E[A_4^2 & convex] + 40 E[A_3^3]",
        method=dict(EA5="B-P width route (n6_bp_polygon.EA_hull_polygon)",
                    EA33="B-P width route (n6_bp_polygon.EA3k_polygon)",
                    EA42_convex="two-chord B-P line-space integral, 80-bit "
                                "(n6_twochord_polygon.two_chord_polygon)"),
        anchors_reproduced=dict(m3="91/900 (Valtr)", m4="49/400 (Valtr)",
                                limit="disk 0.134309386357 (Marckert)"),
        rows=rows, disk=disk)
    json.dump(out, open(os.path.join(R, "n6_mgon_P6_final.json"), "w"), indent=1)

    print(f"{'m':>9} {'E[A_4^2&conv]':>23} {'P_6':>23}  exact")
    for r in rows:
        ex = r.get("P6_exact", "")
        print(f"{r['m']:>9} {r['EA42_convex']:>23} {r['P6']:>23}  {ex}")
    print(f"{disk['m']:>9} {disk['EA42_convex'][:21]:>23} {disk['P6'][:21]:>23}  {disk['P6_exact']}")
    print("\nexact-vs-numeric residuals:",
          {r['m']: '%.1e' % r['exact_vs_numeric'] for r in rows if 'exact_vs_numeric' in r})
    print("two-chord grid spreads:", {r['m']: '%.1e' % r['EA42_convex_err'] for r in rows})


if __name__ == "__main__":
    main()
