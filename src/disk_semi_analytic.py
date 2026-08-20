"""Semi-analytic route (route S) for the DISK, n = 4 and 5.

Ingredients (unit disk, |K| = pi):
 * Renyi-Sulanke edge count. Pairs of points <-> (line, positions): dx1 dx2 = |t1-t2| dt1 dt2 dG,
   integrating positions along a chord of length L (ordered pair) gives int int |t1-t2| = L^3/3. Lines dG = dp dtheta, p in [0,1].
   The chord at distance p has L = 2 sqrt(1-p^2) and cuts off a cap of area
   C(p) = arccos(p) - p sqrt(1-p^2).  Let c = C/pi.
   E[#edges of conv(n pts)] = C(n,2) * (1/pi^2) * 2pi * int_0^1 (L^3/3) [(1-c)^{n-2} + c^{n-2}] dp
   and #edges = #vertices =: N_n for a polygon.
 * Efron: E[A_{n-1}] = 1 - E[N_n]/n   (A = hull area fraction).
 * Buchta identity: P_5 = 1 - 5 E[A_4] + 10 E[A_3^2],  and E[A_3^2] = 3/(32 pi^2) EXACTLY
   (E[det^2] = 3/8 for three uniform points in the unit disk; polynomial moment).
 Hence  P_4 = E[N_4] - 3,   P_5 = E[N_5] - 4 + 15/(16 pi^2).
Anchor: P_4 must equal 1 - 35/(12 pi^2).
"""
from mpmath import mp, mpf, quad, acos, sqrt, pi, binomial, identify, pslq

mp.dps = 40

def EN(n):
    def f(p):
        L3 = (2 * sqrt(1 - p * p)) ** 3
        c = (acos(p) - p * sqrt(1 - p * p)) / pi
        return L3 / 3 * ((1 - c) ** (n - 2) + c ** (n - 2))
    return binomial(n, 2) * (2 * pi / pi ** 2) * quad(f, [0, 1])

EN4 = EN(4); P4 = EN4 - 3
print("E[N_4]  =", EN4)
print("P_4     =", P4)
print("anchor  =", 1 - mpf(35) / (12 * pi ** 2), " diff =", P4 - (1 - mpf(35) / (12 * pi ** 2)))
EN5 = EN(5); P5 = EN5 - 4 + mpf(15) / (16 * pi ** 2)
print("E[N_5]  =", EN5)
print("E[A_4]  =", 1 - EN5 / 5)
print("P_5(disk) =", P5)
# constant recognition attempts (report the basis tried either way)
print("identify(P5, ['pi**2','pi**-2','pi**-4']):", identify(P5, ['pi**2', 'pi**-2', 'pi**-4']))
print("identify(P5, ['pi**-2','pi**-4']):", identify(P5, ['pi**-2', 'pi**-4']))
print("pslq [P5, 1, 1/pi^2, 1/pi^4]:", pslq([P5, 1, 1 / pi ** 2, 1 / pi ** 4], maxcoeff=10 ** 8, maxsteps=10 ** 6))
print("pslq [P5, 1, 1/pi^2, 1/pi^4, 1/pi^6]:", pslq([P5, 1, 1 / pi ** 2, 1 / pi ** 4, 1 / pi ** 6], maxcoeff=10 ** 8, maxsteps=10 ** 6))
