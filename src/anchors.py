"""Exact anchor values for the convex-position probability P_K(n).

Every new code path must reproduce these before any new number from it is
believed (campaign rule, see ../docs/../LESSONS_LEARNED.md in Math/).

Sources (see LITERATURE.md for full citations):
  * Valtr 1995 (DCG 13:637-643): parallelogram  P = [C(2n-2,n-1)/n!]^2
  * Valtr 1996 (Combinatorica 16:567-573): triangle P = 2^n (3n-3)! / [(2n)! ((n-1)!)^3]
  * Sylvester four-point problem, disk/ellipse: P = 1 - 35/(12 pi^2)
  * Sylvester four-point, triangle: 2/3 ; square: 25/36  (special cases of Valtr)
  * 3D five-point (Efron): P = 1 - 5 E[vol tetrahedron]/|K|
      cube: E[V]/|K| = 3977/216000 - pi^2/2160  (Zinani 2003)
      ball: E[V]/|K| = 9/715                    (classical; to be re-verified in LITERATURE.md)
"""
from fractions import Fraction
from math import comb, factorial, pi


def valtr_parallelogram(n: int) -> Fraction:
    """P(n uniform points in a parallelogram are in convex position), Valtr 1995."""
    if n < 3:
        return Fraction(1)
    return Fraction(comb(2 * n - 2, n - 1), factorial(n)) ** 2


def valtr_triangle(n: int) -> Fraction:
    """P(n uniform points in a triangle are in convex position), Valtr 1996."""
    if n < 3:
        return Fraction(1)
    return Fraction(2 ** n * factorial(3 * n - 3), factorial(2 * n) * factorial(n - 1) ** 3)


def sylvester_disk() -> float:
    """P(4 uniform points in a disk/ellipse are in convex position)."""
    return 1 - 35 / (12 * pi ** 2)


def sylvester3d_cube_value() -> float:
    """P(5 uniform points in a cube are in convex position) = 1 - 5*E[V]/|K| (Zinani 2003)."""
    ev = 3977 / 216000 - pi ** 2 / 2160
    return 1 - 5 * ev


def sylvester3d_ball_value() -> float:
    return 1 - 5 * (9 / 715)


ANCHORS = {
    # key: (description, exact-or-high-precision value)
    "square_n4": ("Sylvester, square", float(valtr_parallelogram(4))),          # 25/36
    "square_n5": ("Valtr, square n=5", float(valtr_parallelogram(5))),         # 49/144
    "square_n6": ("Valtr, square n=6", float(valtr_parallelogram(6))),
    "triangle_n4": ("Sylvester, triangle", float(valtr_triangle(4))),          # 2/3
    "triangle_n5": ("Valtr, triangle n=5", float(valtr_triangle(5))),
    "triangle_n6": ("Valtr, triangle n=6", float(valtr_triangle(6))),
    "disk_n4": ("Sylvester, disk", sylvester_disk()),                          # 1 - 35/(12 pi^2)
    "cube3d_n5": ("Sylvester 3D, cube (Zinani)", sylvester3d_cube_value()),
    "ball3d_n5": ("Sylvester 3D, ball (9/715)", sylvester3d_ball_value()),
}

if __name__ == "__main__":
    print("Valtr parallelogram P(n):")
    for n in range(3, 11):
        v = valtr_parallelogram(n)
        print(f"  n={n:2d}: {str(v):>32s} = {float(v):.12f}")
    print("Valtr triangle P(n):")
    for n in range(3, 11):
        v = valtr_triangle(n)
        print(f"  n={n:2d}: {str(v):>32s} = {float(v):.12f}")
    print(f"Sylvester disk n=4: 1 - 35/(12 pi^2) = {sylvester_disk():.12f}")
    print(f"3D cube n=5 (Zinani): {sylvester3d_cube_value():.12f}")
    print(f"3D ball n=5 (9/715): {sylvester3d_ball_value():.12f}")
    assert valtr_parallelogram(4) == Fraction(25, 36)
    assert valtr_parallelogram(5) == Fraction(49, 144)
    assert valtr_triangle(4) == Fraction(2, 3)
    print("anchor sanity asserts OK")


# ---------------------------------------------------------------------------
# Added 2026-08-19 (server session): Marckert 2017 (arXiv:1402.3512) Table (7).
# Exact P_disk(n) for n = 4..8.  Read verbatim from the paper; our Monte Carlo
# agrees with every one of them (|z| <= 1.55, see PROGRESS.md 2026-08-19 00:13).
# ---------------------------------------------------------------------------

def marckert_disk(n: int) -> float:
    """P(n uniform points in a disk are in convex position), Marckert 2017 Table (7)."""
    p = pi
    one_minus = {
        4: 35 / (12 * p ** 2),
        5: 305 / (48 * p ** 2),
        6: (146400 * p ** 2 - 473473) / (11520 * p ** 4),
        7: (512400 * p ** 2 - 2900611) / (23040 * p ** 4),
        8: (62664108221 + 1721664000 * p ** 4 - 18670881600 * p ** 2) / (48384000 * p ** 6),
    }[n]
    return 1 - one_minus


# Alikoski 1939 (via MathWorld "Polygon Triangle Picking"), independently reproduced
# by src/polygon_exact.py (route P) to 50 digits for m = 3..24, 30, 40, 60, 100.
def alikoski_EA3(m: int) -> float:
    """Mean area of a random triangle in a regular m-gon of unit area."""
    from math import cos, sin
    w = 2 * pi / m
    return (9 * cos(w) ** 2 + 52 * cos(w) + 44) / (36 * m ** 2 * sin(w) ** 2)


def p5_regular_mgon(m: int) -> float:
    """P_5(regular m-gon) = 1 - 5(15cos^2 w + 92 cos w + 76)/(36 m^2 sin^2 w), w = 2pi/m.
    (This session: Alikoski + identity (II).  m=3 -> 11/36, m=4 -> 49/144, m->oo -> disk.)"""
    from math import cos, sin
    w = 2 * pi / m
    return 1 - 5 * (15 * cos(w) ** 2 + 92 * cos(w) + 76) / (36 * m ** 2 * sin(w) ** 2)


ANCHORS.update({
    f"disk_n{n}": (f"Marckert 2017, disk n={n}", marckert_disk(n)) for n in (5, 6, 7, 8)
})
