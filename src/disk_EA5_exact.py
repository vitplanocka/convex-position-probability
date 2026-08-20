r"""Exact E[A_5] for the disk (groundwork for n=6), via the c-moment route of THEOREMS.md.

Master formula: E[A_{n-1}] = 1 - (n-1) E[c^{n-2}], and E[N_6] = 15 T_4 = 30 E[c^4], so
    E[A_5](disk) = 1 - 5 E[c^4].
For the unit disk, with p = cos(phi) the cut-off fraction is c = (phi - sin phi cos phi)/pi and
    E[N_n] = C(n,2) * (2/pi) * (8/3) * int_0^{pi/2} sin^4(phi) [(1-c)^{n-2} + c^{n-2}] dphi
(same setup as src/disk_n5_exact.py, whose n=4 case reproduces 1 - 35/(12 pi^2)).
"""
import sympy as sp

phi = sp.symbols('phi', positive=True)
c = (phi - sp.sin(phi) * sp.cos(phi)) / sp.pi


def EN(n):
    integrand = sp.expand(sp.sin(phi) ** 4 * ((1 - c) ** (n - 2) + c ** (n - 2)))
    I = sp.integrate(integrand, (phi, 0, sp.pi / 2))
    return sp.simplify(sp.binomial(n, 2) * (2 / sp.pi) * sp.Rational(8, 3) * I)


for n in (4, 5, 6):
    ENn = EN(n)
    EA = sp.simplify(1 - ENn / n)
    print(f"E[N_{n}](disk) = {ENn}")
    print(f"E[A_{n-1}](disk) = {sp.simplify(EA)}   = {sp.N(EA, 30)}", flush=True)
    if n == 4:
        print("   check P_4 = E[N_4]-3 - (1-35/(12pi^2)) =",
              sp.simplify(ENn - 3 - (1 - sp.Rational(35, 12) / sp.pi ** 2)))
    if n == 6:
        print("   route-P Richardson value (independent): 0.212072074036701")
        print("   difference:", sp.N(EA - sp.Float('0.212072074036701', 20), 20))
