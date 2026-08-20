"""Exact rational verification of route P and identities (I), (II).
Affine invariance lets us reach m = 3, 4, 6 with rational vertices."""
from fractions import Fraction as F
from polygon_exact import exact_rational_moments

CASES = {
    "triangle (0,0),(1,0),(0,1)": [(0, 0), (1, 0), (0, 1)],
    "unit square": [(0, 0), (1, 0), (1, 1), (0, 1)],
    "sheared square (parallelogram)": [(0, 0), (1, 0), (F(4, 3), 1), (F(1, 3), 1)],
    "affinely-regular hexagon": [(2, 0), (1, 3), (-1, 3), (-2, 0), (-1, -3), (1, -3)],
    "random rational pentagon": [(0, 0), (5, F(1, 2)), (7, 4), (3, 6), (-1, 3)],
    "random rational heptagon": [(0, 0), (4, -1), (7, 2), (8, 6), (5, 9), (1, 8), (-2, 4)],
}
EXPECT = {
    "triangle (0,0),(1,0),(0,1)": (F(1, 12), F(2, 3), F(11, 36)),
    "unit square": (F(11, 144), F(25, 36), F(49, 144)),
    "sheared square (parallelogram)": (F(11, 144), F(25, 36), F(49, 144)),
    "affinely-regular hexagon": (F(289, 3888), F(683, 972), F(1373, 3888)),
}

print("=== exact rational route P (Fraction arithmetic, no floating point) ===")
for name, V in CASES.items():
    r = exact_rational_moments(V, nmax=6)
    ea3, ea4, ea5 = r['E[A_3]'], r['E[A_4]'], r['E[A_5]']
    print(f"\n{name}   area={r['area']}")
    print(f"  T_0 = {r['T'][0]}  (must be 2)      E[N_3] = {r['E[N_3]']}  (must be 3)")
    print(f"  E[A_3]   = {ea3} = {float(ea3):.15f}")
    print(f"  E[A_4]   = {ea4}")
    print(f"  IDENTITY (I): E[A_4] - 2 E[A_3] = {ea4 - 2*ea3}   (must be 0)")
    print(f"  E[A_5]   = {ea5} = {float(ea5):.15f}")
    print(f"  P_4      = {r['P_4']} = {float(r['P_4']):.15f}")
    print(f"  E[A_3^2] = {r['E[A_3^2]']} = {float(r['E[A_3^2]']):.15f}")
    print(f"  P_5      = {r['P_5']} = {float(r['P_5']):.15f}")
    if name in EXPECT:
        e3, p4, p5 = EXPECT[name]
        print(f"  CHECK vs literature: E[A_3]-{e3} = {ea3-e3}; P_4-{p4} = {r['P_4']-p4}; P_5-{p5} = {r['P_5']-p5}")
