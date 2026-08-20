"""Numerical check of the fibrewise 'section lemma' proof of the n=5 extremal theorem
(supplied 2026-08-19; see docs/N5_PROOF.md).

Claims tested:
 (L1) For fixed ordered abscissae, R = b1 U1 - b2 U2 + b3 U3 (U_i ~ U[-1,1]) has density
      exactly 1/(2 b2) on [-r, r], r = b2 - b1 - b3 >= 0, and m(t) = E|R+t| satisfies
      m(t) - m(0) = t^2/(2 b2) for |t| <= r.
 (L2) Steiner deficit formula: F(K) - F(SK) = (12/V^5) int_{x1<x2<x3} h1 h2 h3 d^2 (V-b2)/b2 dx,
      tested on K = half-disk {x^2+y^2<=1, y>=0} (vertical sections), SK = ellipse (F = F(disk)).
 (L3) Shaking gap formula: F(ShK) - F(K) = (12/V^5) int h1 h2 h3 (V-b2)/b2 (r^2 - d^2) dx,
      tested on K = disk (d = 0), ShK = half-ellipse ~ half-disk. Must equal the same number.
Reference: F(half-disk) = E[A_3] - E[A_3^2] with E[A_3] = 0.076512497523552 (route P Richardson,
overnight session) and E[A_3^2] = (3/2) det Sigma / V^2 exact; F(disk) = 61/(96 pi^2).
"""
import numpy as np
from math import pi, sqrt
from scipy.stats import qmc

# ---- (L1) plateau + quadratic law on one fibre ------------------------------------------
rng = np.random.default_rng(1)
b1, b2, b3 = 0.3, 1.0, 0.25          # b2 >= b1+b3, r = 0.45
r = b2 - b1 - b3
U = rng.random((3, 4_000_000)) * 2 - 1
R = b1 * U[0] - b2 * U[1] + b3 * U[2]
for t in (-0.4, -0.2, 0.0, 0.2, 0.4):
    dens = np.mean(np.abs(R - t) < 0.01) / 0.02
    m_t = np.mean(np.abs(R + t)); m_0 = np.mean(np.abs(R))
    print(f"(L1) t={t:+.1f}: density {dens:.4f} (claim {1/(2*b2):.4f});  m(t)-m(0) = {m_t-m_0:.5f} (claim {t*t/(2*b2):.5f})")

# ---- reference numbers --------------------------------------------------------------------
V_half = pi / 2
EA3_half = 0.076512497523552
detSigma_half = 0.25 * (0.25 - 16 / (9 * pi ** 2))
EA3sq_half = 1.5 * detSigma_half / V_half ** 2
F_half = EA3_half - EA3sq_half
F_disk = 61 / (96 * pi ** 2)
print(f"E[A_3^2](half-disk) exact = {EA3sq_half:.12f};  F(half-disk) = {F_half:.12f};  P_5(half-disk) = {1-10*F_half:.12f}")
print(f"F(disk) = {F_disk:.12f};  target deficit F(half-disk) - F(disk) = {F_half - F_disk:.12f}")

# ---- (L2) Steiner deficit integral for the half-disk -------------------------------------
def deficit_integral(h, c, V, a=-1.0, b=1.0, m=2**22, seed=7):
    """(12/V^5) int_{a<x1<x2<x3<b} h1 h2 h3 d^2 (V-b2)/b2 dx  via Sobol QMC over the cube + sort."""
    S = qmc.Sobol(d=3, scramble=True, seed=seed).random(m)
    X = np.sort(a + (b - a) * S, axis=1)
    x1, x2, x3 = X[:, 0], X[:, 1], X[:, 2]
    p, q = x2 - x1, x3 - x2
    h1, h2, h3 = h(x1), h(x2), h(x3)
    c1, c2, c3 = c(x1), c(x2), c(x3)
    d = q * c1 - (p + q) * c2 + p * c3
    b2 = (p + q) * h2
    integrand = h1 * h2 * h3 * d * d * (V - b2) / b2
    vol = (b - a) ** 3          # cube volume; ordering handled by sorting (integrand over ordered region = cube/6 * 6 ... )
    # We integrate over the ordered simplex: sorting maps the cube 6-to-1 onto it, so
    # int_simplex f = (1/6) * vol * mean(f(sorted)).
    return 12 / V ** 5 * vol / 6 * integrand.mean(), integrand

hh = lambda x: 0.5 * np.sqrt(np.clip(1 - x * x, 0, None))
cc = lambda x: 0.5 * np.sqrt(np.clip(1 - x * x, 0, None))
vals = [deficit_integral(hh, cc, V_half, seed=s)[0] for s in range(4)]
print(f"(L2) Steiner deficit integral (half-disk -> ellipse): {np.mean(vals):.9f} +- {np.std(vals)/2:.1e}   [target {F_half - F_disk:.9f}]")

# ---- (L3) shaking gap integral for the disk (d = 0) --------------------------------------
def gap_integral_disk(m=2**22, seed=11):
    S = qmc.Sobol(d=3, scramble=True, seed=seed).random(m)
    X = np.sort(-1 + 2 * S, axis=1)
    x1, x2, x3 = X[:, 0], X[:, 1], X[:, 2]
    p, q = x2 - x1, x3 - x2
    h = lambda x: np.sqrt(np.clip(1 - x * x, 0, None))
    h1, h2, h3 = h(x1), h(x2), h(x3)
    b1, b2, b3 = q * h1, (p + q) * h2, p * h3
    r = b2 - b1 - b3
    V = pi
    integrand = h1 * h2 * h3 * (V - b2) / b2 * r * r
    return 12 / V ** 5 * 8 / 6 * integrand.mean(), np.min(r), np.max(b2 / V)
res = [gap_integral_disk(seed=s) for s in range(4)]
print(f"(L3) shaking gap integral (disk -> half-ellipse): {np.mean([x[0] for x in res]):.9f} +- {np.std([x[0] for x in res])/2:.1e}   [target {F_half - F_disk:.9f}]")
print(f"     sanity: min r = {min(x[1] for x in res):.3e} (>= 0 required), max b2/V = {max(x[2] for x in res):.6f} (< 1 required)")
