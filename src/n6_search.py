"""Nelder-Mead min/max search for P_6 over convex-body shape space, plus local
(common-random-numbers) perturbation analysis at the two conjectured extremisers.

The conjecture (open for n >= 6) says
    min_K P_6(K) = P_6(triangle) = 91/900,     max_K P_6(K) = P_6(ellipse) = 0.1343093864.

Two search spaces:
  * k-gons:  vertices q_0 = (0,0), q_1 = (1,0), q_2 = (1/2, sqrt3/2) FIXED (an affine map is
    3-point transitive, so fixing three vertices costs no generality and quotients out the
    6-dimensional affine group), the remaining k-3 vertices free.  Body = conv(all q_i).
  * support functions:  h(theta) = 1 + sum_{k=2..K} (a_k cos k theta + b_k sin k theta),
    which is convex iff h + h'' = 1 - sum k^2 (a_k cos + b_k sin) >= 0.  This is the natural
    space of SMOOTH competitors around the disk (k = 1 is a translation and is excluded).

COMMON RANDOM NUMBERS: mcp6's RNG stream depends only on (seed, nchunk, nsamp), and the
sampler maps a fixed uniform stream continuously into the body, so with a frozen seed the
objective is a nearly-deterministic function of the shape -- nearby shapes share almost all
of their Monte-Carlo error.  Differences of P_6 between nearby bodies are therefore far more
accurate than each value separately (measured: ~40x variance reduction near the disk).
Absolute values quoted anywhere are always re-measured afterwards with fresh seeds.
"""
from __future__ import annotations

import argparse
import json
import math
import time

import numpy as np
from scipy.optimize import minimize

import mcp6
import n6_bodies as NB

P6_TRI, P6_DISK = mcp6.P6_TRI, mcp6.P6_DISK
SQ3 = math.sqrt(3) / 2
BASE = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, SQ3]])


# ------------------------------------------------------------------ parametrisations
def kgon_body(v, k):
    """v: flat array of 2(k-3) free coordinates."""
    P = np.vstack([BASE, np.asarray(v, float).reshape(-1, 2)])
    H = NB.hull(P)
    if len(H) < 3:
        return None
    return H


def support_body(c, K, n=512):
    """c: [a_2,b_2,...,a_K,b_K].  Returns None if h + h'' < 0 anywhere (not convex)."""
    t = np.linspace(0, 2 * math.pi, n, endpoint=False)
    h = np.ones_like(t)
    hpp = np.zeros_like(t)
    hp = np.zeros_like(t)
    for j, k in enumerate(range(2, K + 1)):
        a, b = c[2 * j], c[2 * j + 1]
        h += a * np.cos(k * t) + b * np.sin(k * t)
        hp += k * (-a * np.sin(k * t) + b * np.cos(k * t))
        hpp += -k * k * (a * np.cos(k * t) + b * np.sin(k * t))
    if np.min(h + hpp) <= 1e-9 or np.min(h) <= 1e-9:
        return None
    return np.stack([h * np.cos(t) - hp * np.sin(t), h * np.sin(t) + hp * np.cos(t)], axis=-1)


# ------------------------------------------------------------------ objective
class Obj:
    def __init__(self, build, N, seed, sign=+1, threads=16):
        self.build, self.N, self.seed, self.sign, self.threads = build, N, seed, sign, threads
        self.nev = 0
        self.best = None

    def __call__(self, v):
        V = self.build(v)
        self.nev += 1
        if V is None or len(V) < 3:
            return self.sign * (-1e3)     # penalise degenerate / non-convex
        p = mcp6.mc(V, self.N, self.seed, nthreads=self.threads)["P6"]
        if self.best is None or self.sign * p < self.sign * self.best[0]:
            self.best = (p, np.array(v, float))
        return self.sign * p


def refine(V, N=int(1e9), seeds=(9001, 9002, 9003), threads=16):
    """Re-measure a body with fresh independent seeds."""
    ps6, ps7 = [], []
    for s in seeds:
        r = mcp6.mc(V, N, s, nthreads=threads)
        ps6.append(r["P6"]); ps7.append(r["P7"])
    p6 = float(np.mean(ps6)); p7 = float(np.mean(ps7))
    se6 = math.sqrt(p6 * (1 - p6) / (N * len(seeds)))
    se7 = math.sqrt(p7 * (1 - p7) / (N * len(seeds)))
    return dict(P6=p6, se6=se6, P7=p7, se7=se7, samples=N * len(seeds), seeds=list(seeds),
                spread6=float(np.std(ps6, ddof=1)) if len(ps6) > 1 else 0.0)


# ------------------------------------------------------------------ local perturbation tests
def base_values(V0, N, seeds, threads=16):
    return {s: mcp6.mc(V0, N, s, nthreads=threads)["P6"] for s in seeds}


def crn_delta(base, V1, N, seeds, threads=16):
    """P_6(V1) - P_6(V0) with common random numbers (base = cached P_6(V0) per seed).
    Measured variance reduction vs independent runs: ~14x."""
    d = np.array([mcp6.mc(V1, N, s, nthreads=threads)["P6"] - base[s] for s in seeds])
    return float(d.mean()), float(d.std(ddof=1) / math.sqrt(len(d))) if len(d) > 1 else 0.0


def disk_local(N=int(2e8), seeds=(1, 2, 3, 4, 5, 6), ks=(2, 3, 4, 5, 6, 8), threads=16, nsm=512):
    """Second-order behaviour of P_6 at the disk along h = 1 + eps cos(k theta)."""
    V0 = NB.norm_area(NB.reg(nsm))
    base = base_values(V0, N, seeds, threads)
    rows = []
    for k in ks:
        for eps in (0.25 / k**2, 0.5 / k**2, 0.9 / k**2):
            V1 = NB.norm_area(NB.cos_perturbed_disk(k, eps, nsm))
            m, se = crn_delta(base, V1, N, seeds, threads)
            rows.append(dict(k=k, eps=eps, dP6=m, se=se, dP6_over_eps2=m / eps**2,
                             z=m / se if se else float("nan")))
            print(f"  disk + {eps:.5f} cos({k}t):  dP_6 = {m:+.3e} +- {se:.1e}"
                  f"   (dP_6/eps^2 = {m/eps**2:+.4f})", flush=True)
    return rows


def triangle_local(N=int(2e8), seeds=(1, 2, 3, 4, 5, 6), threads=16):
    """First-order behaviour of P_6 at the triangle: cut one/two/three corners by t, and
    round it by Minkowski-adding a disk of radius r."""
    V0 = NB.norm_area(NB.reg(3))
    base = base_values(V0, N, seeds, threads)
    rows = []
    for kcut in (1, 2, 3):
        for t in (0.002, 0.005, 0.01, 0.02, 0.05):
            V1 = NB.norm_area(NB.truncated_triangle(t, kcut))
            m, se = crn_delta(base, V1, N, seeds, threads)
            rows.append(dict(kind=f"cut{kcut}", t=t, dP6=m, se=se, dP6_over_t2=m / t**2,
                             z=m / se if se else float("nan")))
            print(f"  triangle, {kcut} corner(s) cut by {t:.3f}:  dP_6 = {m:+.3e} +- {se:.1e}"
                  f"   (dP_6/t^2 = {m/t**2:+.4f})", flush=True)
    for r in (0.002, 0.005, 0.01, 0.02, 0.05):
        V1 = NB.norm_area(NB.rounded(NB.reg(3), r))
        m, se = crn_delta(base, V1, N, seeds, threads)
        rows.append(dict(kind="round", t=r, dP6=m, se=se, dP6_over_t2=m / r**2,
                         z=m / se if se else float("nan")))
        print(f"  triangle rounded by r = {r:.3f}:  dP_6 = {m:+.3e} +- {se:.1e}", flush=True)
    return rows


# ------------------------------------------------------------------ drivers
def normalise_to_base(V):
    """Affine map taking V[0],V[1],V[2] to BASE; returns the images of V[3:] (the free
    coordinates of the parametrisation).  Affine maps are 3-point transitive, so every convex
    k-gon has such a representative and the parametrisation is onto shape space."""
    V = np.asarray(V, float)
    M = np.linalg.solve(np.hstack([V[:3], np.ones((3, 1))]), BASE)   # (3x3) -> (3x2)
    return (np.hstack([V[3:], np.ones((len(V) - 3, 1))]) @ M).ravel()


def run_kgon(k, sign, seed, N, nstart, threads, rng):
    out = []
    for st in range(nstart):
        if st == 0:
            v0 = normalise_to_base(NB.reg(k))            # start AT the regular k-gon
        elif st == 1:
            W = NB.hull(np.vstack([NB.reg(k)[:1] * 1.6, NB.reg(k)[1:]]))
            v0 = normalise_to_base(W) if len(W) == k else normalise_to_base(NB.reg(k))
        else:
            v0 = (np.tile([0.5, SQ3 / 3], (k - 3, 1)) + 0.8 * rng.standard_normal((k - 3, 2))).ravel()
        ob = Obj(lambda v, k=k: kgon_body(v, k), N, seed, sign, threads)
        t0 = time.time()
        res = minimize(ob, v0, method="Nelder-Mead",
                       options=dict(maxfev=800, xatol=1e-5, fatol=1e-8, adaptive=True))
        V = kgon_body(res.x, k)
        p, se = (float("nan"), float("nan"))
        if V is not None:
            rr = refine(V, int(5e8), (7001, 7002), threads)
            p, se = rr["P6"], rr["se6"]
        tag = "min" if sign > 0 else "max"
        print(f"  {k}-gon {tag} start {st}: NM {sign*res.fun:.7f} in {ob.nev} evals; "
              f"refined P_6 = {p:.7f} +- {se:.1e}  "
              f"(tri{p-P6_TRI:+.6f}, disk{p-P6_DISK:+.6f})  [{time.time()-t0:.0f}s]", flush=True)
        out.append(dict(k=k, kind=tag, start=st, nm_value=float(sign * res.fun), nev=ob.nev,
                        P6=p, se6=se, nv=int(len(V)) if V is not None else 0,
                        vertices=(NB.norm_area(V).tolist() if V is not None else None)))
    return out


def run_support(K, sign, seed, N, nstart, threads, rng):
    out = []
    for st in range(nstart):
        c0 = 0.15 * rng.standard_normal(2 * (K - 1)) / np.repeat(np.arange(2, K + 1) ** 2, 2)
        ob = Obj(lambda c, K=K: support_body(c, K), N, seed, sign, threads)
        t0 = time.time()
        res = minimize(ob, c0, method="Nelder-Mead",
                       options=dict(maxfev=600, xatol=1e-5, fatol=1e-7, adaptive=True))
        V = support_body(res.x, K)
        p, se = (float("nan"), float("nan"))
        if V is not None:
            rr = refine(NB.norm_area(V), int(5e8), (7001, 7002), threads)
            p, se = rr["P6"], rr["se6"]
        tag = "min" if sign > 0 else "max"
        print(f"  support K={K} {tag} start {st}: NM {sign*res.fun:.7f} in {ob.nev} evals; "
              f"refined P_6 = {p:.7f} +- {se:.1e}  (disk{p-P6_DISK:+.6f})  "
              f"[{time.time()-t0:.0f}s]", flush=True)
        out.append(dict(K=K, kind=tag, start=st, nm_value=float(sign * res.fun), nev=ob.nev,
                        P6=p, se6=se, coeffs=res.x.tolist()))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="all",
                    choices=["all", "local", "kgon", "support"])
    ap.add_argument("--samples", type=float, default=2e7, help="samples per NM evaluation")
    ap.add_argument("--local-samples", type=float, default=2e8)
    ap.add_argument("--seed", type=int, default=4242)
    ap.add_argument("--starts", type=int, default=2)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--out", default="../results/n6_extremal_search.json")
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    N = int(a.samples)
    res = {}
    t00 = time.time()

    if a.mode in ("all", "local"):
        print("=== local perturbation of the DISK (conjectured maximiser) ===", flush=True)
        res["disk_local"] = disk_local(int(a.local_samples), (1, 2, 3, 4, 5, 6), threads=a.threads)
        print("=== local perturbation of the TRIANGLE (conjectured minimiser) ===", flush=True)
        res["triangle_local"] = triangle_local(int(a.local_samples), (1, 2, 3, 4, 5, 6), threads=a.threads)
        json.dump(res, open(a.out, "w"), indent=1)

    if a.mode in ("all", "kgon"):
        print("=== Nelder-Mead over k-gon vertex coordinates ===", flush=True)
        res["kgon"] = []
        for k in (4, 5, 6, 7, 8):
            for sign in (+1, -1):
                res["kgon"] += run_kgon(k, sign, a.seed, N, a.starts, a.threads, rng)
                json.dump(res, open(a.out, "w"), indent=1)

    if a.mode in ("all", "support"):
        print("=== Nelder-Mead over support-function Fourier coefficients ===", flush=True)
        res["support"] = []
        for K in (4, 6, 8):
            for sign in (-1, +1):
                res["support"] += run_support(K, sign, a.seed, N, a.starts, a.threads, rng)
                json.dump(res, open(a.out, "w"), indent=1)

    res["seconds"] = time.time() - t00
    json.dump(res, open(a.out, "w"), indent=1)
    print(f"\ntotal {res['seconds']:.0f}s -> {a.out}")
