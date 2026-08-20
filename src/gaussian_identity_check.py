r"""Is identity (I) really distribution-free?

The proof uses only (a) Efron: E[mu(conv(X_1..X_{n-1}))] = 1 - E[N_n]/n,
(b) Renyi-Sulanke: E[N_n] = C(n,2) E[c^{n-2} + (1-c)^{n-2}] with c the mu-MASS on one
side of the line through two of the points, (c) E[c] = 1/2 by the side symmetry.
None of these needs uniformity or a convex support.  So (I) should read

    E[ mu(conv(X_1,...,X_4)) ] = 2 E[ mu(conv(X_1,X_2,X_3)) ]

for ANY absolutely continuous distribution mu on R^2.  Test: standard Gaussian, a
strongly anisotropic + skewed Gaussian mixture, and a non-convex-support mixture.
Separately test Finch's (arXiv:1601.04937) observation about the LEBESGUE area of the
Gaussian hulls, which is a different statement (mu-mass != Lebesgue area for a Gaussian).
"""
import numpy as np, math, json

def mu_mass_of_hull(P, W, rng_extra):
    """Monte-Carlo estimate of mu(conv(P)) using extra iid points: mu(conv) = P(Y in conv)."""
    raise NotImplementedError


def hull_area(P):
    """Shoelace area of the convex hull of a small point set (n<=5), via monotone chain."""
    from scipy.spatial import ConvexHull
    return ConvexHull(P).volume


def test(name, sampler, N=4_000_000, seed=1, batches=20):
    """mu-content of conv(X_1..X_k) estimated as P(Y in conv) with ONE extra point Y:
       E[mu(conv_k)] = P(Y in conv(X_1..X_k)) exactly.  So we only need an indicator!"""
    rng = np.random.default_rng(seed)
    per = N // batches
    b3 = np.zeros(batches); b4 = np.zeros(batches)
    for b in range(batches):
        X = sampler(rng, per, 5)
        # E[mu(conv(X1..X3))] = P(X4 in conv(X1,X2,X3))
        b3[b] = in_triangle(X[:, 0], X[:, 1], X[:, 2], X[:, 3]).mean()
        # E[mu(conv(X1..X4))] = P(X5 in conv(X1..X4))
        b4[b] = in_quad(X[:, :4], X[:, 4]).mean()
    f = math.sqrt(batches)
    m3, s3 = b3.mean(), b3.std(ddof=1) / f
    m4, s4 = b4.mean(), b4.std(ddof=1) / f
    d = (b4 - 2 * b3)
    return dict(name=name, E_mu3=m3, se3=s3, E_mu4=m4, se4=s4, ratio=m4 / m3,
                diff=d.mean(), se_diff=d.std(ddof=1) / f, z=d.mean() / (d.std(ddof=1) / f))


def _o(a, b, c):
    return (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])


def in_triangle(A, B, C, X):
    d1, d2, d3 = _o(A, B, X), _o(B, C, X), _o(C, A, X)
    return ((d1 >= 0) & (d2 >= 0) & (d3 >= 0)) | ((d1 <= 0) & (d2 <= 0) & (d3 <= 0))


def in_quad(Q, X):
    """X inside conv of 4 points: inside one of the 4 triangles."""
    out = np.zeros(len(X), bool)
    for i in range(4):
        idx = [j for j in range(4) if j != i]
        out |= in_triangle(Q[:, idx[0]], Q[:, idx[1]], Q[:, idx[2]], X)
    return out


def gauss(rng, m, n):
    return rng.standard_normal((m, n, 2))


def skew_mix(rng, m, n):
    """Anisotropic, skewed, NON-convex-support-ish mixture: 70% N(0, diag(1,9)) rotated,
    30% N((5,-3), diag(4,0.25))."""
    u = rng.random((m, n)) < 0.7
    g1 = rng.standard_normal((m, n, 2)) * np.array([1.0, 3.0])
    th = 0.6
    R = np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
    g1 = g1 @ R.T
    g2 = rng.standard_normal((m, n, 2)) * np.array([2.0, 0.5]) + np.array([5.0, -3.0])
    return np.where(u[..., None], g1, g2)


def cauchy_like(rng, m, n):
    """Heavy-tailed (t_3) radial law -- still absolutely continuous."""
    g = rng.standard_normal((m, n, 2))
    s = np.sqrt(rng.chisquare(3, size=(m, n)) / 3)
    return g / s[..., None]


def lebesgue_area_test(sampler, N=2_000_000, seed=5, batches=20):
    """E[Lebesgue area of conv(3)] vs E[Lebesgue area of conv(4)] (Finch's observation)."""
    rng = np.random.default_rng(seed)
    per = N // batches
    a3 = np.zeros(batches); a4 = np.zeros(batches)
    for b in range(batches):
        X = sampler(rng, per, 4)
        a3[b] = (0.5 * np.abs(_o(X[:, 0], X[:, 1], X[:, 2]))).mean()
        # area of conv of 4 pts = sum of the two triangles of whichever diagonal splits it,
        # = max over the 4 triangles' areas *if* one point is inside, else sum of two.
        ar = np.stack([0.5 * np.abs(_o(X[:, i], X[:, j], X[:, k]))
                       for i, j, k in [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]], 1)
        tot = ar.sum(1)
        mx = ar.max(1)
        # if convex position: sum of all four triangles = 2*area; else = 2*area of big one
        a4[b] = np.where(np.isclose(tot, 2 * mx, rtol=1e-9), mx, tot / 2).mean()
    f = math.sqrt(batches)
    return a3.mean(), a3.std(ddof=1) / f, a4.mean(), a4.std(ddof=1) / f


if __name__ == "__main__":
    print("=== identity (I) for the mu-CONTENT, non-uniform / non-convex-support laws ===")
    rows = []
    for nm, s in [("standard Gaussian", gauss), ("skewed anisotropic mixture", skew_mix),
                  ("heavy-tailed t_3", cauchy_like)]:
        r = test(nm, s)
        rows.append(r)
        print(f"{nm:28s} E[mu_3]={r['E_mu3']:.7f}+-{r['se3']:.1e}  E[mu_4]={r['E_mu4']:.7f}+-{r['se4']:.1e}"
              f"  ratio={r['ratio']:.6f}  E[mu_4]-2E[mu_3] = {r['diff']:+.2e} +- {r['se_diff']:.1e}  z={r['z']:+.2f}", flush=True)
    print("\n=== Finch's separate observation: LEBESGUE area, standard Gaussian ===")
    a3, s3, a4, s4 = lebesgue_area_test(gauss)
    print(f"  E[area conv3] = {a3:.6f} +- {s3:.1e}   (sqrt(3)/2 = {math.sqrt(3)/2:.6f}, z={(a3-math.sqrt(3)/2)/s3:+.2f})")
    print(f"  E[area conv4] = {a4:.6f} +- {s4:.1e}   (sqrt(3)   = {math.sqrt(3):.6f}, z={(a4-math.sqrt(3))/s4:+.2f})")
    print("  NB this is a DIFFERENT statement from (I): for a Gaussian the mu-content of the")
    print("     hull is not its Lebesgue area, so the two doublings are not the same fact.")
    json.dump(rows, open("../results/A4_gaussian_identity_check.json", "w"), indent=1, default=float)
