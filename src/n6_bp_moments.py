"""n = 6 groundwork.

(A) Validate the P_6 decomposition and the server's extracted disk moments:
      P_n = sum_{j=0}^{n-3} (-1)^j C(n,j) E[A_{n-j}^j]      (route_moments identity)
    n=6:  P_6 = 1 - 6 E[A_5] + 15 E[A_4^2] - 20 E[A_3^3].
    Disk (Marckert / server extraction):
      E[A_5]   = 7(2400 pi^2 - 3289)/(6912 pi^4)
      E[A_3^3] = 1001/(6400 pi^4)
      E[A_4^2] = (2400 pi^2 + 31031)/(19200 pi^4)
      P_6(disk) = 1 - (146400 pi^2 - 473473)/(11520 pi^4)   [Marckert 2017 Table (7)]

(B) An EXACT route to the triangle-area moments E[A_3^k] for ANY convex body, via the
    two-point Blaschke-Petkantschin (base-line) formula.  For P1,P2,P3 iid uniform in K,
    Delta = area(P1 P2 P3) = (1/2)|t1 - t2| * d(P3, line(P1,P2)), and the planar 2-point
    B-P identity dP1 dP2 = |t1 - t2| dt1 dt2 dp dphi gives

      \int_{K^3} Delta^k = (2^{1-k}/((k+2)(k+3))) \int_0^pi \int_p w(p,phi)^{k+3}
                            [ \int_t |t - p|^k w(t,phi) dt ] dp dphi,

    where w(.,phi) is the width function of K in normal direction phi (the length of the
    slice of K by the line {x . n_phi = .}).  Then E[A_3^k] = (1/V^{k+3}) \int_{K^3} Delta^k.
    This is route P's philosophy (integrate over line space) applied to an area MOMENT
    rather than an edge count; k=2 is polynomial, but k=3 (an absolute third moment) is the
    new object needed for n=6.  Here it is validated on the disk against the known
    k=1,2,3 values; the polygon version is in n6_bp_polygon.py.
"""
import mpmath as mp

mp.mp.dps = 40


# ---------------------------------------------------------------- (A) framework
def check_framework():
    import sympy as sp
    pi = sp.pi
    EA5 = 7 * (2400 * pi**2 - 3289) / (6912 * pi**4)
    EA33 = sp.Rational(1001) / (6400 * pi**4)
    EA42 = (2400 * pi**2 + 31031) / (19200 * pi**4)
    P6_marckert = 1 - (146400 * pi**2 - 473473) / (11520 * pi**4)
    P6_assembled = 1 - 6 * EA5 + 15 * EA42 - 20 * EA33
    diff = sp.simplify(P6_assembled - P6_marckert)
    print("[A] P_6 decomposition vs Marckert disk value")
    print("    assembled  P_6(disk) =", sp.nsimplify(sp.simplify(P6_assembled)))
    print("    Marckert   P_6(disk) =", sp.simplify(P6_marckert))
    print("    difference           =", diff, " ->", "PASS" if diff == 0 else "FAIL")
    # also the derived convex-quadrilateral second moment
    EA42_conv = sp.simplify(EA42 - 4 * EA33)
    print("    E[A_4^2 & convex](disk) = E[A_4^2] - 4 E[A_3^3] =", EA42_conv, "=", sp.nsimplify(EA42_conv))
    return diff == 0


# ---------------------------------------------------------------- (B) disk B-P
def disk_width(t):
    if abs(t) >= 1:
        return mp.mpf(0)
    return 2 * mp.sqrt(1 - t * t)


def EA3k_disk_bp(k):
    """E[A_3^k] for the unit disk via the B-P width-function integral (phi trivial)."""
    V = mp.pi
    coef = mp.mpf(2) ** (1 - k) / ((k + 2) * (k + 3))

    def slab(p):  # \int_t |t-p|^k w(t) dt
        return mp.quad(lambda t: abs(t - p) ** k * disk_width(t), [-1, p, 1])

    def outer(p):
        return disk_width(p) ** (k + 3) * slab(p)

    integral_pphi = mp.pi * mp.quad(outer, [-1, 0, 1])   # \int_0^pi dphi = pi (disk symmetric)
    IK3 = coef * integral_pphi                            # = \int_{K^3} Delta^k
    return IK3 / V ** (k + 3)


def check_disk_bp():
    known = {
        1: mp.mpf(35) / (48 * mp.pi**2),
        2: mp.mpf(3) / (32 * mp.pi**2),
        3: mp.mpf(1001) / (6400 * mp.pi**4),
    }
    print("[B] disk E[A_3^k] via Blaschke-Petkantschin width integral")
    ok = True
    for k in (1, 2, 3):
        got = EA3k_disk_bp(k)
        exp = known[k]
        err = abs(got - exp)
        ok &= err < mp.mpf(10) ** (-30)
        print(f"    k={k}: BP = {mp.nstr(got, 25)}   known = {mp.nstr(exp, 25)}   |err| = {mp.nstr(err, 3)}")
    print("    ->", "PASS" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    a = check_framework()
    print()
    b = check_disk_bp()
    print()
    print("SUMMARY:", "framework", "OK" if a else "FAIL", "| disk B-P route", "OK" if b else "FAIL")
