r"""A3: independent Monte-Carlo confirmation of identities (I) and (II) on nine bodies.

For each body:
  * DIRECT tester (convex_position.py, numba monotone-chain hull; tester A cross-check
    on a subsample)         -> P_4^MC, P_5^MC
  * hull-area route (route_moments.hull_areas, Jarvis march, no shared code)
                            -> E[A_3]^MC, E[A_4]^MC, E[A_3^2]^MC
  * route P / exact         -> E[A_3]^P (exact for polygons, Richardson for smooth bodies),
                               E[A_3^2] exact from det Sigma
Tests
  (I)  E[A_4] - 2 E[A_3] = 0                              -> z_I
  (II) P_5 = 1 - 10 E[A_3] + 10 E[A_3^2]                   -> z_II  (route-P prediction vs direct MC)
       P_4 = 1 - 4 E[A_3]                                  -> z_IV
"""
import json, math, sys, time
import numpy as np

import convex_position as CP
import route_moments as RM
from polygon_exact import moments_fast, regular_polygon, exact_rational_moments
from bodies_extra import SAMPLERS, AREAS, POLY, exact_EA3sq
from route_p_bodies import body_route_p
from anchors import valtr_parallelogram, valtr_triangle

# register the extra bodies with both engines
for nm, s in SAMPLERS.items():
    CP.BODIES[nm] = (s, 2)
    RM.BODY_AREA[nm] = AREAS[nm]

POLY_VERTS = {
    "square": [(0, 0), (1, 0), (1, 1), (0, 1)],
    "triangle": [(0, 0), (1, 0), (0, 1)],
    "pentagon": regular_polygon(5), "hexagon": regular_polygon(6), "octagon": regular_polygon(8),
}
BODIES = ["square", "triangle", "disk", "pentagon", "hexagon", "octagon",
          "ellipse3", "halfdisk", "stadium"]


def reference(body):
    """(E[A_3], E[A_3^2]) from route P / exact, plus a label for how it was obtained."""
    if body in ("square", "triangle"):
        r = exact_rational_moments(POLY_VERTS[body], nmax=5)
        return float(r["E[A_3]"]), float(r["E[A_3^2]"]), "exact rational"
    if body in ("pentagon", "hexagon", "octagon"):
        r = moments_fast(POLY_VERTS[body], nmax=5)
        return r["E[A_3]"], r["E[A_3^2]"], "route P exact polygon"
    r = body_route_p(body, verbose=False, nmax=5)
    return r["E[A_3]"], r["E[A_3^2]"], "route P Richardson"


def mc_moments(body, samples, seed, batches=25):
    sampler, _ = CP.BODIES[body]
    area = RM.BODY_AREA[body]
    rng = np.random.default_rng(seed)
    per = samples // batches
    b3 = np.zeros((batches, 2))
    b4 = np.zeros(batches)
    for b in range(batches):
        s3 = np.zeros(2); s4 = 0.0; done = 0
        while done < per:
            m = min(2_000_000, per - done)
            A3 = RM.hull_areas(sampler(rng, m, 3)) / area
            A4 = RM.hull_areas(sampler(rng, m, 4)) / area
            s3 += [A3.sum(), (A3 ** 2).sum()]
            s4 += A4.sum()
            done += m
        b3[b] = s3 / per
        b4[b] = s4 / per
    f = math.sqrt(batches)
    return {"E[A_3]": b3[:, 0].mean(), "se_E[A_3]": b3[:, 0].std(ddof=1) / f,
            "E[A_3^2]": b3[:, 1].mean(), "se_E[A_3^2]": b3[:, 1].std(ddof=1) / f,
            "E[A_4]": b4.mean(), "se_E[A_4]": b4.std(ddof=1) / f,
            "E[A_4]-2E[A_3]": (b4 - 2 * b3[:, 0]).mean(),
            "se_diff": (b4 - 2 * b3[:, 0]).std(ddof=1) / f}


if __name__ == "__main__":
    NP = int(float(sys.argv[1])) if len(sys.argv) > 1 else 200_000_000   # direct-tester samples
    NM = int(float(sys.argv[2])) if len(sys.argv) > 2 else 100_000_000   # moment samples
    out = []
    for body in BODIES:
        t0 = time.time()
        ea3, ea3sq, how = reference(body)
        mm = mc_moments(body, NM, seed=101)
        r4 = CP.estimate(body, 4, NP, seed=202, both=False)
        r5 = CP.estimate(body, 5, NP, seed=303, both=False)
        r5x = CP.estimate(body, 5, 4_000_000, seed=404, both=True)   # A/B cross-check
        P4_pred = 1 - 4 * ea3
        P5_pred = 1 - 10 * ea3 + 10 * ea3sq
        zI = mm["E[A_4]-2E[A_3]"] / mm["se_diff"]
        zIV = (r4["p_hat"] - P4_pred) / r4["std_err"]
        zII = (r5["p_hat"] - P5_pred) / r5["std_err"]
        zA3 = (mm["E[A_3]"] - ea3) / mm["se_E[A_3]"]
        zA3sq = (mm["E[A_3^2]"] - ea3sq) / mm["se_E[A_3^2]"]
        row = dict(body=body, ref_method=how, EA3_ref=ea3, EA3sq_ref=ea3sq,
                   P4_pred=P4_pred, P5_pred=P5_pred,
                   P4_mc=r4["p_hat"], se_P4=r4["std_err"], P5_mc=r5["p_hat"], se_P5=r5["std_err"],
                   samples_direct=NP, samples_moments=NM,
                   testers_agree=r5x.get("testers_agree"),
                   z_identityI=zI, z_identityII=zII, z_P4=zIV, z_EA3=zA3, z_EA3sq=zA3sq,
                   seconds=time.time() - t0, **mm)
        out.append(row)
        print(f"{body:9s} [{how}]  E[A_3]={ea3:.12f}", flush=True)
        print(f"          (I)   E[A_4]-2E[A_3] = {mm['E[A_4]-2E[A_3]']:+.3e} +- {mm['se_diff']:.1e}   z={zI:+.2f}")
        print(f"          (P4)  pred {P4_pred:.9f}  MC {r4['p_hat']:.9f} +- {r4['std_err']:.1e}  z={zIV:+.2f}")
        print(f"          (II)  pred {P5_pred:.9f}  MC {r5['p_hat']:.9f} +- {r5['std_err']:.1e}  z={zII:+.2f}"
              f"   testersAgree={r5x.get('testers_agree')}  ({time.time()-t0:.0f}s)", flush=True)
    json.dump(out, open("../results/A3_identity_mc.json", "w"), indent=1, default=float)
    zs = [max(abs(r["z_identityI"]), abs(r["z_identityII"]), abs(r["z_P4"])) for r in out]
    print(f"\nmax |z| over all bodies and all three tests: {max(zs):.2f}")
    print("wrote ../results/A3_identity_mc.json")
