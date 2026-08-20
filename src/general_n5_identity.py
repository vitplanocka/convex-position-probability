"""General planar identities (derived 2026-08-18; to be verified numerically here):

  (I)  E[A_4] = 2 E[A_3]   for every planar convex body K
       (Renyi-Sulanke edge counts E[N_n] = C(n,2) E[(1-c)^{n-2} + c^{n-2}], c = area
        fraction on the left of the directed chord through two random points, E[c] = 1/2
        by side-symmetry  =>  E[N_4] = 12 E[c^2],  E[N_5] = 30 E[c^2] - 5 = (5/2) E[N_4] - 5;
        Efron: E[A_{n-1}] = 1 - E[N_n]/n   =>  E[A_4] = 2 E[A_3]).
  (II) P_5(K) = 1 - 5 E[A_4] + 10 E[A_3^2]                 (Buchta-type identity)
             = 1 - 10 E[A_3] + 10 E[A_3^2] = 1 - 10 E[A_3 (1 - A_3)]
             = (5/2) P_4(K) - 3/2 + 10 E[A_3^2],
       and E[A_3^2] = (3/2) det(Sigma_K) / |K|^2  (E[det^2] = 6 det Sigma for three iid
       points with covariance Sigma; A_3 = |det|/(2|K|)).
Exact checks: triangle 11/36 and square 49/144 (Valtr), disk 1 - 305/(48 pi^2).
This script: MC of E[A_3], E[A_4], E[A_3^2] for square/triangle/disk/pentagon/hexagon,
plus the covariance-determinant value of E[A_3^2] for polygons.
"""
import math
import numpy as np
from convex_position import BODIES
from route_moments import hull_areas, BODY_AREA
from anchors import valtr_parallelogram, valtr_triangle


def polygon_cov_det(V):
    """(det Sigma, area) for the uniform distribution on the polygon with vertices V (ccw),
    via fan triangulation from the origin (origin must be inside or use signed areas)."""
    V = np.asarray(V, float)
    n = len(V)
    A_tot = 0.0
    M1 = np.zeros(2)
    M2 = np.zeros((2, 2))
    o = np.zeros(2)
    for i in range(n):
        a, b = V[i], V[(i + 1) % n]
        T = 0.5 * (a[0] * b[1] - a[1] * b[0])  # signed
        g = (a + b + o) / 3
        m2 = (np.outer(a, a) + np.outer(b, b) + np.outer(o, o) + np.outer(a + b + o, a + b + o)) / 12
        A_tot += T
        M1 += T * g
        M2 += T * m2
    mean = M1 / A_tot
    S = M2 / A_tot - np.outer(mean, mean)
    return np.linalg.det(S), A_tot


def regular_polygon(k):
    ang = 2 * np.pi * np.arange(k) / k
    return np.stack([np.cos(ang), np.sin(ang)], -1)


def exact_EA3sq(body):
    if body == 'disk':
        return 3 / (32 * math.pi ** 2)
    if body == 'square':
        return 1 / 96
    if body == 'triangle':
        return 1 / 72
    if body in ('pentagon', 'hexagon', 'octagon'):
        k = {'pentagon': 5, 'hexagon': 6, 'octagon': 8}[body]
        d, A = polygon_cov_det(regular_polygon(k))
        return 1.5 * d / A ** 2
    raise KeyError(body)


def mc(body, samples=2_000_000, seed=5, batch=200_000):
    sampler, _ = BODIES[body]
    rng = np.random.default_rng(seed)
    area = BODY_AREA[body]
    A3 = np.concatenate([hull_areas(sampler(rng, batch, 3)) / area for _ in range(samples // batch)])
    A4 = np.concatenate([hull_areas(sampler(rng, batch, 4)) / area for _ in range(samples // batch)])
    n3, n4 = len(A3), len(A4)
    return (A3.mean(), A3.std() / math.sqrt(n3), (A3 ** 2).mean(), (A3 ** 2).std() / math.sqrt(n3),
            A4.mean(), A4.std() / math.sqrt(n4))


if __name__ == "__main__":
    print("exact E[A_3^2] (covariance-determinant route):")
    for b in ['square', 'triangle', 'disk', 'pentagon', 'hexagon']:
        print(f"  {b:9s} {exact_EA3sq(b):.9f}")
    print("square should be 1/96 =", 1 / 96, ", triangle 1/72 =", 1 / 72, ", disk 3/(32pi^2) =", 3 / (32 * math.pi ** 2))
    P4_exact = {'square': 25 / 36, 'triangle': 2 / 3, 'disk': 1 - 35 / (12 * math.pi ** 2)}
    P5_ref = {'square': float(valtr_parallelogram(5)), 'triangle': float(valtr_triangle(5)),
              'disk': 1 - 305 / (48 * math.pi ** 2)}
    for body in ['square', 'triangle', 'disk', 'pentagon', 'hexagon']:
        m3, s3, m3sq, s3sq, m4, s4 = mc(body)
        print(f"{body:9s} E[A3]={m3:.6f}+-{s3:.1e}  E[A4]={m4:.6f}+-{s4:.1e}  E[A4]/E[A3]={m4 / m3:.5f}  "
              f"E[A3^2]={m3sq:.6f}+-{s3sq:.1e} (exact {exact_EA3sq(body):.6f})")
        P4 = P4_exact.get(body, 1 - 4 * m3)
        P5 = 2.5 * P4 - 1.5 + 10 * exact_EA3sq(body)
        ref = P5_ref.get(body)
        print(f"          P5 via identity = {P5:.7f}" + (f"   exact ref {ref:.7f}" if ref else "   (P4 from MC here)"))
