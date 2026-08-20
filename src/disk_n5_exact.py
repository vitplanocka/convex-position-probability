"""Exact symbolic evaluation of P_5(disk) via the semi-analytic route.
E[N_n] = C(n,2) * (2/pi) * int_0^1 (8/3)(1-p^2)^{3/2} [(1-c)^{n-2} + c^{n-2}] dp,
c(p) = (arccos p - p sqrt(1-p^2))/pi.   Substitute p = cos(phi):
E[N_n] = C(n,2) * (2/pi)*(8/3) * int_0^{pi/2} sin^4(phi) [(1-c)^{n-2} + c^{n-2}] dphi,
c = (phi - sin(phi)cos(phi))/pi.
P_4 = E[N_4] - 3 ;  P_5 = E[N_5] - 4 + 15/(16 pi^2)  (E[A_3^2] = 3/(32 pi^2) exact).
"""
import sympy as sp
phi = sp.symbols('phi', positive=True)
c = (phi - sp.sin(phi) * sp.cos(phi)) / sp.pi
def EN(n):
    integrand = sp.expand(sp.sin(phi) ** 4 * ((1 - c) ** (n - 2) + c ** (n - 2)))
    I = sp.integrate(integrand, (phi, 0, sp.pi / 2))
    return sp.simplify(sp.binomial(n, 2) * (2 / sp.pi) * sp.Rational(8, 3) * I)
EN4 = EN(4); P4 = sp.simplify(EN4 - 3)
print("E[N_4] =", EN4, "   P_4 =", P4, "   check vs 1-35/(12pi^2):", sp.simplify(P4 - (1 - sp.Rational(35, 12) / sp.pi ** 2)))
EN5 = EN(5); P5 = sp.simplify(EN5 - 4 + sp.Rational(15, 16) / sp.pi ** 2)
print("E[N_5] =", EN5)
print("E[A_4](disk) =", sp.simplify(1 - EN5 / 5))
print("P_5(disk) =", P5, "  =", sp.N(P5, 30))
print("check vs 1-305/(48pi^2):", sp.simplify(P5 - (1 - sp.Rational(305, 48) / sp.pi ** 2)))
# also E[det^2]=3/8 check for E[A_3^2]
x1,y1,x2,y2,x3,y3,r,t = sp.symbols('x1 y1 x2 y2 x3 y3 r t', real=True)
det = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
# moments of uniform unit disk: E[x^2]=1/4, E[x]=E[xy]=0; independence -> E[det^2] = 3*E[x^2]*E[(y2-y3)^2] = 3*(1/4)*(1/2)
Ex2 = sp.integrate(sp.integrate((r*sp.cos(t))**2 * r, (r,0,1)), (t,0,2*sp.pi)) / sp.pi
print("E[x^2] over unit disk =", Ex2, " => E[det^2] =", 3*Ex2*(2*Ex2), " E[A_3^2] = E[det^2]/4/pi^2 =", sp.simplify(3*Ex2*(2*Ex2)/4/sp.pi**2))
