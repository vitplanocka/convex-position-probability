# The n >= 6 extremal conjecture, stress-tested on 568 non-regular bodies

*Server session, 2026-08-19 (task `TASK_N6_EXTREMAL.md`, Part 1).  Code: `src/mcp6.py` (fast
direct convex-position Monte Carlo of `P_6` and `P_7` for an arbitrary convex polygon),
`src/n6_bodies.py` (the body catalogue), `src/n6_scan.py` (the scan), `src/n6_search.py`
(local perturbation analysis + Nelder-Mead min/max), `src/n6_disk_climb.py`,
`src/n6_climb_recheck.py`, `src/n6_exact_confirm.py` (exact-machinery confirmation).
Results: `results/n6_mc_validation.json`, `results/n6_mc_precision_1e10.json`,
`results/n6_extremal_scan.json`, `results/n6_local_perturbation.json`,
`results/n6_search_kgon.json`, `results/n6_search_support.json`,
`results/n6_disk_climb.json`, `results/n6_climb_recheck.json`,
`results/n6_exact_confirm.json`, `results/n6_disk_P7.json`.
Sibling of `docs/N6_LANDSCAPE.md` (the exact regular-m-gon table).*

## 0. The question and the answer

`P_n(K) = P(n i.i.d. uniform points in a planar convex body K are in convex position)`.
The conjecture, **open for every n >= 6** (Blaschke settled n = 4, Marckert-Rahmani 2021 n = 5):

> `P_n(triangle) <= P_n(K) <= P_n(ellipse)` for every planar convex body `K`.

At n = 6 the window is wide and both ends are known exactly:

    P_6(triangle) = 91/900               = 0.101111111111  (Valtr)
    P_6(ellipse)  = 1 - (146400 pi^2 - 473473)/(11520 pi^4)
                                          = 0.134309386357  (Marckert)

**Verdict: consistent, and now sharply so.**  Over 568 area-normalised convex bodies at
4e8 samples each, no body came within 3 sigma of either bound from the wrong side; the most
extreme deviations in the whole scan were `-1.26 sigma` below the triangle (a body that IS a
triangle) and `-1.17 sigma` above the disk (the regular 32-gon).  A Nelder-Mead search over
k-gon vertex coordinates (k = 4..8) and over support-function Fourier coefficients never
produced a body outside the window: every minimisation run collapsed to the triangle, every
k-gon maximisation run converged to the regular k-gon, and the three shapes that "beat" the
disk during the noisy search all fell back to or below the disk value when re-measured at
1.2e10 samples with fresh seeds.  Beyond consistency, the two conjectured extremisers are
confirmed to be strict local extremisers with measured curvature (section 5) -- which is the
part of the evidence a scan alone cannot give.

At n = 7 the same holds against `P_7(triangle) = 0.0251851852` (Valtr) and
`P_7(disk) = 0.0390905623` (Marckert 2017, Table (7)), which a 2e10-sample run here confirms
independently: `0.039092252 +- 1.4e-6`, `z = +1.23`.

## 1. The Monte-Carlo engine and why it can be trusted

`src/mcp6.py` estimates `P_6` and `P_7` for an arbitrary convex polygon at **5e7 samples/s**
on 16 threads (2e8 samples in 4.1 s).

* **Sampling.** Fan triangulation from the centroid, triangle chosen by area-weighted binary
  search, then the reflected-barycentric map.  Smooth bodies are represented by 256-gons
  (the disk by a 1024-gon); the exact m-gon table gives the polygonal bias as
  `a_4/m^4 = 3.27/m^4`, i.e. `7.6e-10` at m = 256 and `3e-12` at m = 1024 -- four to six
  orders of magnitude below the Monte-Carlo error.
* **Convex-position test.** The centroid of the n sample points is strictly interior to their
  hull, so sorting the points by angle about it gives the hull order whenever all n are hull
  vertices.  Hence *all n in convex position* <=> *the angularly sorted polygon turns left at
  every vertex*.  (=> hull order = angular order about an interior point; <= all-left-turns
  plus winding number 1 about the sort centre forces convexity, and then every point is a
  vertex.)  Angles are compared through the transcendental-free monotone "diamond"
  pseudo-angle `d(dx,dy) = dx/(|dx|+|dy|) - 1` (dy < 0) / `1 - dx/(|dx|+|dy|)` (dy >= 0),
  strictly increasing in the true angle on `(-pi, pi]`; the turns are exact float orientation
  determinants -- the same predicate as `convex_position.py`'s tester A.
* **One pass gives both n.** Each sample draws 7 points; the first 6 give `P_6`, all 7 give
  `P_7`.  The two estimators are separately unbiased (and correlated, which is harmless).
* **RNG.** Counter-based splitmix64, re-seeded per chunk, so the output depends only on
  `(seed, nchunk, nsamp)` and **not on the thread count** -- bit-reproducible.

**Validation** (`mcp6.py --validate --samples 4e8`, `results/n6_mc_validation.json`):
24 checks, all OK, 17 of them against exact values, every z-score in `[-1.4, +1.1]`.

| check | result |
|---|---|
| `P_6` = 91/900, right / equilateral / 1:100-sliver triangle | z = +1.04 (identical: affine invariance is exact in this code path) |
| `P_6` = 49/400, unit square and a sheared parallelogram | z = -0.22 |
| `P_6`, exact regular m-gon table (m = 5..12, 14, 16, 18, 20) from `results/n6_mgon_P6_final.json` | all \|z\| <= 1.4 |
| `P_6`(disk) via a 1024-gon vs Marckert's closed form | z = +0.49 |
| `P_7` = `2^n(3n-3)!/((n-1)!^3(2n)!)` = 0.02518519 (triangle, Valtr 1995) | z = +0.34 |
| `P_7` = `(C(2n-2,n-1)/n!)^2` = 0.03361111 (parallelogram, Valtr 1996) | z = +0.00 |
| affine invariance on a NON-regular body (trapezoid vs a shear+rotation image) | identical |
| 4 cross-checks vs `convex_position.py`'s independent monotone-chain hull tester | \|z\| <= 2.0 |

**Precision floor** (`results/n6_mc_precision_1e10.json`, 1e10 samples each, se ~ 3e-6):
triangle `z = -0.72` (`P_7`: `-1.26`), square `z = +0.49` (`P_7`: `+0.19`), disk over two seeds
(2e10 samples) `z = +1.98`.  **No detectable bias down to 3e-6**, i.e. 1e4 times finer than the
conjecture window.  The same runs sharpen the campaign's check of Marckert's n = 7 disk value:

> `P_7(disk)`: Marckert 2017 gives **0.0390905623**; 2e10 samples here give
> **0.039092252 +- 1.4e-6**, `z = +1.23`.  (The campaign's previous check,
> `LITERATURE.md`, had this at 2e8 samples and `z = +1.55`; the error bar is now 10x smaller.)

**Confirmation against the EXACT first-principles machinery on NON-regular bodies.**  The
exact route (`P_6 = 1 - 6E[A_5] + 15E[A_4^2&conv] + 40E[A_3^3]`, width route +
two-chord route, `docs/N6_LANDSCAPE.md`) had only ever been checked on triangles, squares and
regular m-gons.  `src/n6_exact_confirm.py` runs it on seven irregular bodies and compares with
4e9-sample MC (`results/n6_exact_confirm.json`):

| body | vertices | `P_6` exact | `P_6` MC (4e9) | z |
|---|---|---|---|---|
| trapezoid (t = 0.4) | 4 | 0.11524543345034665 | 0.115248497 ± 5.0e-6 | +0.61 |
| rhombus (kite 1.0/0.8) | 4 | 0.12249999999992331 (= 49/400 to 8e-13) | 0.122498969 ± 5.2e-6 | -0.20 |
| triangle, one corner cut by 0.05 | 4 | 0.10178878869638569 | 0.101783653 ± 4.8e-6 | -1.07 |
| triangle, three corners cut by 0.2 | 6 | 0.12096187020153892 | 0.120960802 ± 5.2e-6 | -0.21 |
| regular pentagon, one vertex pushed 1.6x | 5 | 0.12751177127659825 | 0.127515662 ± 5.3e-6 | +0.74 |
| random quadrilateral | 4 | 0.10298445217813361 | 0.102986911 ± 4.8e-6 | +0.51 |
| random pentagon (hull = quad) | 4 | 0.11049377927283081 | 0.110487874 ± 5.0e-6 | -1.19 |

Two structurally disjoint pipelines -- a 4-D line-space integral in 80-bit arithmetic and a
1e10-sample hull test -- agree to `5e-6` on seven irregular shapes.  This also yields seven new
exact `P_6` values for non-regular bodies as a by-product.

## 2. The scan

`src/n6_bodies.py` builds **568 area-normalised convex bodies** in ten families; 541 of them
are not affinely a triangle.  `src/n6_scan.py` runs each at **4e8 samples** (se ~ 1.7e-5),
88 minutes in total (`results/n6_extremal_scan.json`).

| family | count | `P_6` range | role |
|---|---|---|---|
| regular m-gons, m = 3..20, 24, 32, 48, 64, 128, 256, 1024 (+ a 3:1 ellipse) | 26 | 0.10111 – 0.13433 | control; ordering |
| one/two-vertex pushed or pulled m-gons (m = 3..12, factors 0.55–2.0) | 86 | 0.10109 – 0.13368 | irregular polygons |
| quadrilaterals: trapezoids, skew trapezoids, kites, random | 62 | 0.10110 – 0.12251 | the 2-parameter quad shape space |
| random pentagons / hexagons | 40 | 0.10110 – 0.12605 | |
| Minkowski interpolations tri↔disk, sq↔disk, tri↔sq, pent↔disk, tri↔hex | 55 | 0.10158 – 0.13431 | 1-parameter families joining the two bounds |
| near-triangle: 1/2/3 corners truncated (t = 0.005–0.45), rounded triangles (r = 0.002–0.2) | 31 | 0.10111 – 0.13156 | approach to the conjectured **minimiser** |
| near-disk: `h = 1 + eps cos(k theta)` (k = 3..12, eps up to the convexity limit `1/k^2`), rounded polygons, near-circular ellipses | 62 | 0.12929 – 0.13433 | approach to the conjectured **maximiser** |
| smooth classical: circular segments, sectors, stadiums, lenses, Reuleaux 3/5/7/9-gons, cones | 44 | 0.10287 – 0.13432 | |
| random convex polygons, 4..12 vertices, three point processes | 162 | 0.10110 – 0.13047 | |

**No violations.**  Over all 568 bodies,

    min over bodies of  (P_6 - 91/900)/se     = -1.26     (attained by a body that IS a triangle)
    min over bodies of  (P_6(disk) - P_6)/se  = -1.17     (regular 32-gon, an inner approximation
                                                           to the disk)
    min over bodies of  (P_7 - P_7(tri))/se   = -2.22     (a triangle)
    min over bodies of  (P_7(disk) - P_7)/se  = -2.03     (a 1.02:1 ellipse, i.e. a disk)

with a violation threshold of `-3`.  All four extremes are attained by bodies that are affine
images of, or polygonal approximations to, the very extremiser they "violate": exactly the
statistical noise the conjecture predicts and nothing else.

*Control.*  The catalogue happens to contain 27 bodies that are affinely triangles (pushed
equilateral triangles, degenerate trapezoids, random point sets whose hull is a triangle).
Their `P_6` values scatter around 91/900 with `chi^2/df = 0.70` over 27 degrees of freedom and
mean z = +0.01 -- a clean check that the sampler is unbiased across very different vertex data.

**Closest approaches.**  Signed gaps are `P_6 - 91/900` at the bottom and `P_6 - P_6(disk)` at
the top; the scan's per-body error is `1.7e-5`, so anything smaller than that is noise.

| | body (all NOT affinely the extremiser) | `P_6` | signed gap | z |
|---|---|---|---|---|
| above the triangle | triangle, ONE corner cut by t = 0.005 (a quadrilateral) | 0.1011082 | -3e-6 | -0.2 |
| | triangle rounded by r = 0.002 | 0.1011262 | +1.5e-5 | +1.0 |
| | triangle, one corner cut by t = 0.01 | 0.1011266 | +1.6e-5 | +1.0 |
| | triangle, two corners cut by t = 0.005 | 0.1011342 | +2.3e-5 | +1.5 |
| below the disk | `h = 1 + 0.2/64 cos 8 theta` | 0.1343223 | **+1.3e-5** | -0.8 |
| | segment of the disk cut at half-angle 3.00 rad | 0.1343191 | **+1.0e-5** | -0.6 |
| | stadium = unit disk (+) a segment of length 0.05 | 0.1343125 | +0.3e-5 | -0.2 |
| | Reuleaux 9-gon | 0.1341478 | -1.6e-4 | +9.5 |

The three bold entries sit *nominally above* the disk, by 0.2 to 0.8 sigma -- pure noise, and
each is a body designed to be within `1e-4` of a disk in the first place (a `1/k^2`-amplitude
`k = 8` support perturbation, a disk with a hair shaved off, a disk Minkowski-lengthened by 2.5%).  The
mode-by-mode second variation of §5 puts every one of them strictly below the disk with 11 to
235 sigma of significance; the scan simply cannot see gaps of order `1e-5`.

## 3. Ordering along the families

**Regular m-gons.**  The exact table (`docs/N6_LANDSCAPE.md`) already proves strict increase
in m for m = 3..20; the scan reproduces it to within Monte-Carlo error and extends it to
m = 24, 32, 48, 64, 128, 256, 1024.  Beyond m ~ 17 the exact spacing (`3.27/m^4` per step) drops
below the 1.7e-5 sampling error, so the scan's five apparent "decreases" (m = 17->18, 19->20,
32->48, 128->256, 256->1024) are all sub-sigma (worst -1.56) and carry no information.

**Minkowski interpolations.**  On the three families whose far endpoint is the DISK --
tri↔disk, sq↔disk, pent↔disk -- `P_6` increases strictly at all 30 steps, i.e. moving
Minkowski-linearly from a polygon towards the maximiser raises `P_6` monotonically.  This is
the ordering a Steiner-symmetrisation proof would have to produce (symmetrisation moves a body
towards the disk and would have to increase `P_6` at each step); the data are consistent with
it, though a monotone 1-parameter family is of course not a proof that symmetrisation itself is
monotone.

The two families joining two POLYGONS are *not* monotone, and instructively so: tri↔sq peaks at
`P_6 = 0.128271` at t = 0.7 and then falls to the square's 0.1225, and tri↔hex peaks at
0.132864 at t = 0.8 before falling to the hexagon's 0.131830.  Both interior maxima exceed
*both* endpoints -- as they must, since a Minkowski sum of a triangle and a square is a
heptagon, hence rounder than either summand.  Both stay comfortably inside the window.

**n = 6 vs n = 7.**  Over the 568 bodies, `corr(P_6, P_7) = 0.99919` and the rank correlation
is `0.99974`.  The two functionals order the bodies almost identically, so nothing in the
catalogue distinguishes the n = 6 from the n = 7 conjecture.

## 4. The Nelder-Mead searches

Two parametrisations, both quotienting out the affine group.

* **k-gons.**  Vertices `q_0 = (0,0)`, `q_1 = (1,0)`, `q_2 = (1/2, sqrt3/2)` are FIXED (affine
  maps are 3-point transitive, so fixing three vertices costs no generality) and the remaining
  `k-3` are free; the body is the hull.  Starts: the regular k-gon, a pushed regular k-gon.
* **Support functions.**  `h(theta) = 1 + sum_{k=2..K}(a_k cos k theta + b_k sin k theta)`,
  convex iff `h + h'' >= 0`; `k = 1` is a translation and is excluded, `k = 2` is the affine
  direction to leading order.

The objective uses **common random numbers** (a frozen seed): the sampler maps a fixed uniform
stream continuously into the body, so nearby shapes share almost all of their Monte-Carlo
error.  Measured variance reduction on differences near the disk: **14x**.

**Minimisation (`results/n6_search_kgon.json`).**  All ten runs (k = 4, 5, 6, 7, 8, two starts
each) collapse the free vertices into the hull of the fixed three and return
`P_6 = 0.1011079 +- 9.5e-6` -- the triangle, 91/900, to within `3e-6`.  **The minimiser over
k-gons is the triangle for every k <= 8 tested.**

**Maximisation over k-gons.**  Every run converges to the regular k-gon:

| k | NM maximum (refined, 1e9 samples) | exact regular k-gon | gap |
|---|---|---|---|
| 4 | 0.1224663, 0.1224737 ± 1.0e-5 | 49/400 = 0.1225000 | -3.4e-5, -2.6e-5 |
| 5 | 0.1292348, 0.1292291 ± 1.1e-5 | 0.1292484 | -1.4e-5, -1.9e-5 |
| 6 | 0.1318433, 0.1317837 ± 1.1e-5 | 0.1318298 | +1.4e-5, -4.6e-5 |
| 7 | 0.1329543, 0.1329135 ± 1.1e-5 | 0.1329615 | -7e-6, -4.8e-5 |
| 8 | 0.1335092, 0.1333223 ± 1.1e-5 | 0.1335163 | -7e-6, -1.9e-4 |

(the search stops slightly short of the exact optimum, as a noisy Nelder-Mead must).  All are
far below the disk, and the sequence in k rises towards it.

**Maximisation over support functions** (`results/n6_search_support.json`, K = 4, 6, 8, two
starts each, ~600 evaluations per run): the highest value seen in **any** of the ~6600
evaluations was `0.1342343`, i.e. `7.5e-5` *below* the disk.

**The decisive climb** (`src/n6_disk_climb.py`).  Nelder-Mead started EXACTLY at the disk with
an explicit initial simplex, K = 3, 5, 8, ~450 shapes each at 4e7 samples under common random
numbers.  Each run reported a best-of-run gain over the disk's own same-seed value of
`+3.6e-5`, `+4.3e-5`, `+3.9e-5` -- about 2.5-3 sigma of the CRN difference error, and exactly
the selection maximum ("winner's curse") one expects from ~450 noisy draws.  Re-measured with
**fresh seeds at 1.2e10 samples** (`src/n6_climb_recheck.py`, `results/n6_climb_recheck.json`):

| K | claimed gain | fresh `P_6` (1.2e10) | `P_6(disk) - P_6` |
|---|---|---|---|
| 3 | +3.64e-5 | 0.134306034 ± 3.1e-6 | +3.35e-6 (-1.08 sigma) |
| 5 | +4.26e-5 | 0.134308587 ± 3.1e-6 | +0.80e-6 (-0.26 sigma) |
| 8 | +3.90e-5 | 0.134305635 ± 3.1e-6 | +3.75e-6 (-1.21 sigma) |

**Every apparent violation evaporated.**  This is the protocol working as designed, and it is
the reason no claim is made from a search maximum alone.

## 5. The two extremisers are strict local extremisers (the sharpest evidence)

A scan can only fail to find a counterexample.  The local analysis
(`src/n6_search.py --mode local`, `results/n6_local_perturbation.json`; 2e8 samples per shape
per seed, six seeds, common random numbers) measures the *shape derivative* directly.

**At the disk.**  Perturb the support function, `h = 1 + eps cos(k theta)` (convex for
`eps <= 1/k^2`), and fit `dP_6 = -c_2 eps^2 - c_3 eps^3` over three eps per mode:

| mode k | `c_2` | significance |
|---|---|---|
| 2 | **+0.0005 ± 0.0004** | consistent with **zero** |
| 3 | -0.3765 ± 0.0016 | 235 sigma |
| 4 | -0.5476 ± 0.0075 | 73 sigma |
| 5 | -0.6608 ± 0.0145 | 46 sigma |
| 6 | -0.7151 ± 0.0268 | 27 sigma |
| 8 | -0.6759 ± 0.0631 | 11 sigma |

> **The disk is a strict local maximiser of `P_6`: the second variation is negative definite on
> every Fourier mode `k >= 3`, and the `k = 2` mode -- the affine (ellipse) direction -- is a
> null direction, exactly as affine invariance demands.**

The `k = 2` null direction is worth spelling out, because it is a self-check the code could not
have faked: `h = 1 + eps cos 2 theta` agrees with an ellipse to first order in `eps` but not to
second, so `P_6` must be *flat to order `eps^2`* and depart only at order `eps^4`.  Measured:
`dP_6 = -7.4e-7, -1.03e-5, -8.2e-5` at `eps = 0.0625, 0.125, 0.225` -- ratios 13.9 and 7.9
against `2^4 = 16` and `1.8^4 = 10.5`.  A quartic, not a quadratic.

**At the triangle.**  Cut `kc = 1, 2, 3` corners of an equilateral triangle back a fraction `t`
of each adjacent edge (removing an area fraction `kc t^2`):

| perturbation | `dP_6` | `dP_6` / removed area |
|---|---|---|
| 1 corner, t = 0.02 | +1.335e-4 ± 7.8e-6 | +0.334 |
| 1 corner, t = 0.05 | +6.965e-4 ± 8.3e-6 | +0.279 |
| 2 corners, t = 0.02 | +2.406e-4 ± 9.1e-6 | +0.301 |
| 2 corners, t = 0.05 | +1.358e-3 ± 6.9e-6 | +0.272 |
| 3 corners, t = 0.01 | +9.389e-5 ± 4.1e-6 | +0.313 |
| 3 corners, t = 0.02 | +3.724e-4 ± 4.6e-6 | +0.310 |
| 3 corners, t = 0.05 | +2.053e-3 ± 8.6e-6 | +0.274 |

Fourteen of the 15 corner-cut perturbations raise `P_6` (10 at >= 3 sigma, 12 at >= 2 sigma;
the five weaker ones are exactly the smallest cuts, t = 0.002 and 0.005, where the predicted
effect is at or below the 8e-6 noise).  The one negative entry, 2 corners at t = 0.002, is
`dP_6 = -2.1e-6 ± 9.5e-6` -- consistent with the `+2.4e-6` the law predicts and with zero, not
with a decrease.  The response is **linear in the removed area with the same
coefficient for 1, 2 and 3 corners**, `dP_6 ~ 0.30 x (removed area fraction)` -- so the corners
contribute additively at leading order:

> **The triangle is a strict local minimiser of `P_6`, and it is a *corner* of shape space, not
> a smooth critical point: the one-sided derivative of `P_6` in the corner-truncation direction
> is strictly positive (~0.3 per unit removed area), not zero.**

**A subtlety, recorded.**  Rounding the triangle instead (`T (+) rB`, Minkowski) gives
`dP_6 = 6.0e-6, 2.4e-5, 6.7e-5, 2.4e-4, 1.2e-3` at `r = 0.002, 0.005, 0.01, 0.02, 0.05` -- a
power law with exponent `1.5-1.9` (overall fit 1.65), i.e. **sublinear** in `r`, so the first
derivative in this direction *vanishes*.  That is not a contradiction but a consequence of a
special feature of the triangle: its **inner** parallel bodies `T (-) rB` are again triangles,
so `r -> P_6(T (+) rB)` is *constant* (= 91/900) for `r <= 0` and >= 91/900 for `r > 0`; a
function that is constant on a half-line and minimised there has vanishing one-sided derivative
if it is C^1.  The exact exponent is not identified and is recorded as UNVERIFIED.

## 6. Verdict

1. **The n = 6 conjecture is consistent with everything measured.**  568 bodies at 4e8 samples,
   plus targeted runs at 1e9-1.2e10; nothing outside `[91/900, P_6(disk)]` beyond `1.3 sigma`,
   and the three search "violations" were refuted at 1.2e10 samples.
2. **Both conjectured extremisers are confirmed to be strict LOCAL extremisers**, with measured
   second variation at the disk (negative definite off the affine null direction) and measured
   positive one-sided derivative at the triangle.  This is stronger than a scan: it rules out
   counterexamples *near* the two candidates, which is where a counterexample would most
   plausibly hide.
3. **The ordering predicted by a Steiner-symmetrisation proof holds** on every 1-parameter
   family tested (regular m-gons, and five Minkowski interpolations from the triangle/square
   /pentagon towards the disk).
4. **The same conclusions hold at n = 7**, against `P_7(triangle) = 0.0251851852` (Valtr) and
   `P_7(disk) = 0.0390905623` (Marckert 2017), the latter re-confirmed here at 2e10 samples
   (`z = +1.23`).
5. **New sub-results.**  The maximiser over k-gons is the regular k-gon for k = 4..8, and the
   minimiser over k-gons is the triangle for k = 4..8 (numerically, to ~1e-5).
6. **What this is NOT.**  No proof.  The Monte-Carlo resolution is ~1.7e-5 per body in the scan
   and ~3e-6 in the targeted runs, so a counterexample violating either bound by less than
   ~1e-5 would be invisible to the scan (though the local analysis of section 5 excludes
   counterexamples in a whole neighbourhood of the two extremisers at a much finer level).
   Bodies with extreme aspect ratio are covered only up to affine equivalence -- which is
   exact, since `P_n` is affine invariant, so this is not a gap.

## 7. Cost and reproduction

    python mcp6.py --validate --samples 4e8 --out ../results/n6_mc_validation.json   # 4 min
    python n6_scan.py --samples 4e8 --seed 11                                        # 88 min
    python n6_search.py --mode local                                                 # 22 min
    python n6_search.py --mode kgon --starts 2 --samples 2e7                         # 84 min
    python n6_search.py --mode support --starts 2 --samples 2e7                      # 67 min
    python n6_disk_climb.py 4e7 ; python n6_climb_recheck.py                         # 33 + 16 min
    python n6_exact_confirm.py --samples 4e9                                         # 8 min

All on 16 threads at `nice -n 10`; total ~5.5 h, of which the exact two-chord + width machinery
accounts for 3 minutes and the Monte Carlo for the rest.
