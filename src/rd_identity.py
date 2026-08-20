r"""R^d analogues of  E[A_4] = 2 E[A_3]  (d = 2)  and  E[A_5] = (5/2) E[A_4]  (d = 3).

MAIN CLAIM (proved here by a Radon-circuit argument, checked deterministically and by MC):

    For any convex body K in R^d and P_1,...,P_{d+2} i.i.d. uniform in K,

        vol conv(P_1,...,P_{d+2})  =  (1/2) sum_{i=1}^{d+2} vol conv(P_1,..,^P_i,..,P_{d+2})

    POINTWISE (almost surely), hence in expectation

        E[A_{d+2}] = ((d+2)/2) E[A_{d+1}]        for EVERY d >= 1.

Proof.  d+2 points in general position in R^d form a circuit: the affine dependence
sum_i lambda_i P_i = 0, sum_i lambda_i = 0 is unique up to scale and has all lambda_i != 0.
Let P = {i : lambda_i > 0}, N = {i : lambda_i < 0} be the Radon partition.  A circuit has
exactly two triangulations, T_P = { conv(all \ {i}) : i in P } and T_N = { conv(all \ {i}) :
i in N }, each covering conv(all) with disjoint interiors.  Hence
sum_{i in P} D_i = sum_{i in N} D_i = vol conv(all), where D_i = vol conv(all \ {i}), and
adding the two gives sum_{i=1}^{d+2} D_i = 2 vol conv(all).  []

Consequences / relation to the Dehn-Sommerville route.  The route used in this campaign for
d = 2, 3 (Efron  E[f_0(n)] = n(1 - E[A_{n-1}])  +  the facet formula
E[f_{d-1}(n)] = C(n,d) E[(1-c)^{n-d} + c^{n-d}]  +  a linear relation f_{d-1} = alpha f_0 + beta)
STOPS at d = 3: the f-vectors of simplicial d-polytopes have floor(d/2) degrees of freedom
modulo Dehn-Sommerville, which is 1 exactly for d = 2, 3, so f_0 and f_{d-1} are
Dehn-Sommerville-independent for every d >= 4 (in d = 4, f_2 = 2 f_3 and f_1 = f_0 + f_3, and
f_1 is not a c-moment).  The identity nevertheless survives in all dimensions -- by the circuit
argument, not by Dehn-Sommerville.  So the R^3 note's conjecture "the phenomenon survives in ODD
dimensions" is wrong in both directions: the DS derivation dies at d = 4 (even AND odd), while
the identity itself is dimension-free.

Checks in this file:
  --circuit   pointwise identity by Qhull volumes, d = 2..7, plus the Radon-partition
              triangulation statement (both triangulations reproduce the hull volume)
  --ds        f-vector experiment in R^4: f_0 and f_3 are not functionally related
  --mc        E[A_{d+2}] = ((d+2)/2) E[A_{d+1}] by Monte Carlo in R^3 and R^4 on several bodies,
              with the two sides computed by INDEPENDENT code (simplex determinants vs Qhull)
"""
from __future__ import annotations

import argparse
import json
import math
import time
from itertools import combinations

import numpy as np
from scipy.spatial import ConvexHull


def simplex_vol(P):
    """P: (..., d+1, d) -> volume of the simplex."""
    d = P.shape[-1]
    M = P[..., 1:, :] - P[..., :1, :]
    return np.abs(np.linalg.det(M)) / math.factorial(d)


def radon(P):
    """P: (d+2, d). Returns lambda with sum lambda = 0, sum lambda_i P_i = 0."""
    d = P.shape[1]
    A = np.vstack([P.T, np.ones(d + 2)])          # (d+1) x (d+2)
    u, s, vt = np.linalg.svd(A)
    return vt[-1]


def check_circuit(d, ntrial=400, rng=None, tol=1e-9):
    rng = rng or np.random.default_rng(0)
    worst_id = 0.0
    worst_tri = 0.0
    bad = 0
    for _ in range(ntrial):
        P = rng.standard_normal((d + 2, d))
        try:
            hv = ConvexHull(P).volume
        except Exception:
            continue
        D = np.array([simplex_vol(np.delete(P, i, axis=0)) for i in range(d + 2)])
        worst_id = max(worst_id, abs(0.5 * D.sum() - hv) / hv)
        lam = radon(P)
        Pp = np.where(lam > 0)[0]
        Nn = np.where(lam < 0)[0]
        if len(Pp) == 0 or len(Nn) == 0:
            bad += 1
            continue
        worst_tri = max(worst_tri, abs(D[Pp].sum() - hv) / hv, abs(D[Nn].sum() - hv) / hv)
    return dict(d=d, ntrial=ntrial, worst_rel_identity=worst_id,
                worst_rel_radon_triangulation=worst_tri, degenerate=bad,
                ok=bool(worst_id < tol and worst_tri < tol))


def check_ds(d, ntrial=3000, rng=None):
    """Exhibit random simplicial d-polytopes with the SAME f_0 but DIFFERENT f_{d-1}, which
    proves f_{d-1} is not a function of f_0 (hence no Dehn-Sommerville relation between them).
    Predicted: possible exactly when floor(d/2) >= 2, i.e. for every d >= 4."""
    rng = rng or np.random.default_rng(1)
    seen = {}
    pairs = []
    for _ in range(ntrial):
        k = int(rng.integers(d + 2, d + 9))
        P = rng.standard_normal((k, d))
        try:
            h = ConvexHull(P)
        except Exception:
            continue
        seen.setdefault(len(h.vertices), set()).add(len(h.simplices))
    for f0, sset in sorted(seen.items()):
        if len(sset) > 1:
            pairs.append((f0, sorted(sset)))
    return dict(d=d, dof_mod_dehn_sommerville=d // 2,
                f0_to_ftop={k: sorted(v) for k, v in sorted(seen.items())},
                counterexamples=pairs,
                ftop_is_function_of_f0=not bool(pairs))


# --------------------------------------------------------------------- bodies for the MC check
def bodies(d):
    """dict name -> sampler(rng, m) -> (m, d) uniform points."""
    def ball(rng, m):
        g = rng.standard_normal((m, d))
        g /= np.linalg.norm(g, axis=1, keepdims=True)
        return g * rng.random((m, 1)) ** (1.0 / d)

    def cube(rng, m):
        return rng.random((m, d))

    def simplex(rng, m):
        e = -np.log(rng.random((m, d + 1)))
        e /= e.sum(axis=1, keepdims=True)
        return e[:, :d]

    def crosspolytope(rng, m):
        # uniform in the l1 ball: Dirichlet on the simplex with random signs
        e = -np.log(rng.random((m, d + 1)))
        e /= e.sum(axis=1, keepdims=True)
        x = e[:, :d]
        return x * rng.choice([-1.0, 1.0], size=(m, d))

    def cylinder(rng, m):
        g = rng.standard_normal((m, d - 1))
        g /= np.linalg.norm(g, axis=1, keepdims=True)
        r = rng.random((m, 1)) ** (1.0 / (d - 1))
        return np.hstack([g * r, rng.random((m, 1))])

    def cone(rng, m):
        h = rng.random((m, 1)) ** (1.0 / d)
        g = rng.standard_normal((m, d - 1))
        g /= np.linalg.norm(g, axis=1, keepdims=True)
        r = rng.random((m, 1)) ** (1.0 / (d - 1))
        return np.hstack([g * r * (1 - h), h])

    return dict(ball=ball, cube=cube, simplex=simplex, crosspoly=crosspolytope,
                cylinder=cylinder, cone=cone)


def mc_identity(d, name, sampler, nsamp=200_000, seed=3, batch=20_000):
    """E[A_{d+1}] by simplex determinants (numpy) and E[A_{d+2}] by Qhull; also the pointwise
    circuit identity on every sample."""
    rng = np.random.default_rng(seed)
    # volume of K by MC of a known reference?  Not needed: A_k is a RATIO, and we get |K| from
    # the same sampler using the (d+2)-point identity only in ratio form.  So compute |K| exactly.
    s1 = 0.0
    s2 = 0.0
    q1 = 0.0
    q2 = 0.0
    n = 0
    worst = 0.0
    keep = np.array([[j for j in range(d + 2) if j != i] for i in range(d + 2)])
    while n < nsamp:
        m = min(batch, nsamp - n)
        P = sampler(rng, m * (d + 2)).reshape(m, d + 2, d)
        v1 = simplex_vol(P[:, : d + 1, :])                    # the (d+1)-point simplex
        Dall = simplex_vol(P[:, keep, :])                      # (m, d+2) leave-one-out volumes
        hv = np.array([ConvexHull(P[s]).volume for s in range(m)])
        worst = max(worst, float(np.max(np.abs(0.5 * Dall.sum(axis=1) - hv) / hv)))
        s1 += float(v1.sum()); q1 += float((v1 * v1).sum())
        s2 += float(hv.sum()); q2 += float((hv * hv).sum())
        n += m
    m1, m2 = s1 / n, s2 / n
    se1 = math.sqrt(max(q1 / n - m1 * m1, 0) / n)
    se2 = math.sqrt(max(q2 / n - m2 * m2, 0) / n)
    ratio = m2 / m1
    return dict(d=d, body=name, nsamp=n, EV_dp1=m1, EV_dp2=m2, se_dp1=se1, se_dp2=se2,
                ratio=ratio, se_ratio=ratio * math.hypot(se1 / m1, se2 / m2),
                predicted=(d + 2) / 2.0, worst_rel_pointwise=worst)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--circuit", action="store_true")
    ap.add_argument("--ds", action="store_true")
    ap.add_argument("--mc", action="store_true")
    ap.add_argument("--nsamp", type=float, default=1e5)
    ap.add_argument("--out", default="../results/rd_identity.json")
    a = ap.parse_args()
    if not (a.circuit or a.ds or a.mc):
        a.circuit = a.ds = a.mc = True
    out = {}
    rng = np.random.default_rng(20260819)

    if a.circuit:
        print("=== pointwise circuit identity  vol conv = (1/2) sum_i D_i ===")
        out["circuit"] = []
        for d in range(2, 8):
            r = check_circuit(d, 400, rng)
            out["circuit"].append(r)
            print(f"  d={d}: worst rel |(1/2)sum D_i - vol| = {r['worst_rel_identity']:.2e};"
                  f"  Radon triangulations {r['worst_rel_radon_triangulation']:.2e}"
                  f"   {'OK' if r['ok'] else 'FAIL'}")

    if a.ds:
        print("\n=== is f_{d-1} a function of f_0? (Dehn-Sommerville independence) ===")
        out["ds"] = []
        for d in (2, 3, 4, 5, 6):
            r = check_ds(d, 3000, rng)
            out["ds"].append(r)
            ex = r["counterexamples"][:1]
            print(f"  d={d}: floor(d/2) = {d//2};  f_{d-1} determined by f_0: "
                  f"{r['ftop_is_function_of_f0']}"
                  + (f"   e.g. f_0 = {ex[0][0]} occurs with f_{d-1} in {ex[0][1]}" if ex else ""))

    if a.mc:
        print("\n=== MC: E[A_(d+2)] / E[A_(d+1)] vs (d+2)/2 ===")
        out["mc"] = []
        for d in (3, 4, 5):
            for name, s in bodies(d).items():
                t0 = time.time()
                r = mc_identity(d, name, s, int(a.nsamp), seed=17 + d)
                r["seconds"] = round(time.time() - t0, 1)
                out["mc"].append(r)
                print(f"  d={d} {name:10s} ratio = {r['ratio']:.6f} +- {r['se_ratio']:.1e}"
                      f"  predicted {r['predicted']:.4f}"
                      f"   z = {(r['ratio']-r['predicted'])/r['se_ratio']:+.2f}"
                      f"   pointwise worst {r['worst_rel_pointwise']:.1e}  [{r['seconds']}s]",
                      flush=True)
    json.dump(out, open(a.out, "w"), indent=1, default=float)
    print(f"\n-> {a.out}")
