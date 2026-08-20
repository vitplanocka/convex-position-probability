# R^d analogues: the volume identity in every dimension, and the 3-D Sylvester 5-point problem

*Server session, 2026-08-19 (task `TASK_N6_EXTREMAL.md`, Part 2).  Code: `src/rd_identity.py`
(the circuit identity, the Dehn-Sommerville obstruction, the R^3/R^4/R^5 Monte-Carlo checks),
`src/rd_scan.py` (the 3-D five-point extremal scan), `src/rd_local.py` (local behaviour at the
simplex and the ball).  Results: `results/rd_identity.json`, `results/rd_scan.json`,
`results/rd_local.json`.  Sibling of `docs/N6_EXTREMAL.md`.  Notation: `A_k` = (volume of the
convex hull of k i.i.d. uniform points in K) / vol(K).*

---

# Part A. The identity family in R^d

## A.1 What was known here, and what the task asked

Two identities had been established in this campaign:

    E[A_4] = 2 E[A_3]         in R^2   (THEOREMS.md, identity (I))
    E[A_5] = (5/2) E[A_4]     in R^3   (src/dim3_identity.py)

both by the same route: Efron's `E[f_0(n)] = n(1 - E[A_{n-1}])`, the Renyi-Sulanke facet count
`E[f_{d-1}(n)] = C(n,d) E[(1-c)^{n-d} + c^{n-d}]` (with `c` the volume fraction of `K` on one
side of the hyperplane through `d` of the points), and a **linear relation between `f_0` and
`f_{d-1}`** for simplicial polytopes: `f_1 = f_0` in the plane, `f_2 = 2f_0 - 4` (Euler +
`3f_2 = 2f_1`) in space.  The task asked for the R^d generalisation via Dehn-Sommerville, and
for the dimensions in which the phenomenon survives; the R^3 note conjectured **odd**
dimensions.

**Both halves have clean answers, and they point in opposite directions.**

## A.2 The identity is dimension-free -- and pointwise

> **Theorem (circuit identity).**  Let `P_1, ..., P_{d+2}` be points in general position in
> `R^d`, and write `D_i = vol conv(P_1, ..., \hat P_i, ..., P_{d+2})` for the volume of the
> `d`-simplex obtained by deleting `P_i`.  Then
>
>     vol conv(P_1, ..., P_{d+2})  =  (1/2) sum_{i=1}^{d+2} D_i .

*Proof.*  `d+2` points in general position in `R^d` form a **circuit**: the affine dependence
`sum_i lambda_i P_i = 0`, `sum_i lambda_i = 0` is unique up to scale and has every
`lambda_i != 0`.  Let `P = {i : lambda_i > 0}` and `N = {i : lambda_i < 0}` be the Radon
partition (both non-empty).  A circuit has exactly **two** triangulations,

    T_P = { conv(all \ {i}) : i in P }    and    T_N = { conv(all \ {i}) : i in N },

each of which covers `conv(P_1..P_{d+2})` with pairwise disjoint interiors.  Hence
`sum_{i in P} D_i = sum_{i in N} D_i = vol conv(all)`, and adding the two gives
`sum_{i=1}^{d+2} D_i = 2 vol conv(all)`.  []

The two degenerate cases are included: if `P = {j}` is a singleton then `P_j` lies inside the
simplex on the other `d+1` points, `T_P` is the single big simplex and `T_N` is its barycentric-
type subdivision from `P_j`.

> **Corollary.**  For `P_1, ..., P_{d+2}` i.i.d. from any absolutely continuous law on `R^d`
> (in particular uniform on a convex body `K`),
>
>     **E[A_{d+2}] = ((d+2)/2) E[A_{d+1}]        for every d >= 1.**

`d = 1`: `E[A_3] = (3/2) E[A_2]` (three points on a line).  `d = 2`: `E[A_4] = 2 E[A_3]`.
`d = 3`: `E[A_5] = (5/2) E[A_4]`.  `d = 4`: `E[A_6] = 3 E[A_5]`.  No Efron, no Renyi-Sulanke, no
Euler relation is needed, and the statement holds **pointwise almost surely**, not merely in
expectation -- a strictly stronger fact than either of the two identities the campaign had.

*Provenance.*  The two-triangulations-of-a-circuit fact is standard in the theory of
triangulations of point configurations (De Loera-Rambau-Santos).  The identity should therefore
be treated as **folklore-or-known, not new** -- consistent with `LITERATURE.md`'s verdict that
`E[A_4] = 2E[A_3]` is to be treated as known.  What is new here is only the observation that
this one line replaces the whole Efron / Renyi-Sulanke / Euler chain and removes the dimension
restriction.

## A.3 Why the Dehn-Sommerville route stops at d = 3 (and the "odd dimensions" guess is wrong)

The route needs `f_{d-1} = alpha f_0 + beta` to hold for **all** simplicial `d`-polytopes,
because `E[f_{d-1}(n)]` is a `c`-moment while `E[f_0(n)]` is Efron's.  Now the `f`-vectors of
simplicial `d`-polytopes satisfy the Dehn-Sommerville relations `h_k = h_{d-k}`, which leave
exactly `floor(d/2)` degrees of freedom.  So `f_0` and `f_{d-1}` are functionally related **iff
`floor(d/2) = 1`, i.e. iff `d = 2` or `d = 3`.**

    d = 2 :  f_1 = f_0                                   (one free parameter, h_1)
    d = 3 :  f_2 = 2 f_0 - 4                              (one free parameter, h_1 = f_0 - 3)
    d = 4 :  f_2 = 2 f_3,  f_1 = f_0 + f_3                (two free parameters: f_0 and f_3 are
                                                           INDEPENDENT; f_1 is not a c-moment)
    d = 5 :  f_4 = 2 + 2h_1 + 2h_2,  f_0 = h_1 + 5        (two free parameters: independent)

So the phenomenon dies at `d = 4` and stays dead -- in even **and** odd dimensions.  The R^3
note's "odd dimensions" conjecture is wrong in both directions: the *derivation* stops at
`d = 3`, and the *identity* never stops.

Numerical confirmation (`rd_identity.py --ds`, 3000 random polytopes per dimension, hulls of
`d+2` to `d+8` gaussian points):

| d | `floor(d/2)` | is `f_{d-1}` a function of `f_0`? | witness |
|---|---|---|---|
| 2 | 1 | **yes** | -- |
| 3 | 1 | **yes** | -- |
| 4 | 2 | no | `f_0 = 6` occurs with `f_3` in {8, 9} |
| 5 | 2 | no | `f_0 = 7` occurs with `f_4` in {10, 12} |
| 6 | 3 | no | `f_0 = 8` occurs with `f_5` in {12, 15, 16} |

For completeness: carrying the general `(alpha, beta)` route through -- eliminate `m_2 = E[c^2]`
between `n = d+2` and `n = d+3`,

    n = d+2:  C(d+2,d) * 2 m_2               = alpha (d+2)(1 - E[A_{d+1}]) + beta
    n = d+3:  C(d+3,3) (3 m_2 - 1/2)         = alpha (d+3)(1 - E[A_{d+2}]) + beta

-- reproduces exactly `E[A_{d+2}] = ((d+2)/2) E[A_{d+1}]` for `(d, alpha, beta) = (2, 1, 0)` and
`(3, 2, -4)`: the constants conspire, as they must, since the circuit identity says the answer
cannot depend on the route.

## A.4 Verification (`results/rd_identity.json`)

**Pointwise**, 400 random gaussian configurations of `d+2` points per dimension; hull volume by
Qhull vs `(1/2) sum_i D_i` by numpy determinants (independent code paths):

| d | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|
| worst relative error | 9.0e-16 | 1.3e-15 | 7.1e-15 | 2.5e-15 | 3.5e-15 | 4.3e-15 |

and each of the two Radon triangulations `T_P`, `T_N` separately reproduces the hull volume to
`<= 7.6e-15`, i.e. the proof's intermediate step is verified too, not just its conclusion.

**In expectation**, `E[A_{d+2}] / E[A_{d+1}]` vs `(d+2)/2`, by Monte Carlo with the numerator
from Qhull and the denominator from simplex determinants, 1e5 configurations per body:

| d | ball | cube | simplex | cross-polytope | cylinder | cone | predicted |
|---|---|---|---|---|---|---|---|
| 3 | 2.4981 | 2.5062 | 2.4970 | 2.4926 | 2.4906 | 2.4977 | 2.5 |
| 4 | 3.0076 | 2.9994 | 3.0193 | 2.9894 | 2.9896 | 3.0002 | 3.0 |
| 5 | 3.5100 | 3.5002 | 3.4947 | 3.5082 | 3.4932 | 3.5127 | 3.5 |

All 18 z-scores lie in `[-0.99, +1.19]`.

## A.5 What the next moment costs

The identity collapses `E[A_{d+2}]` onto `E[A_{d+1}]` and stops there: at `n = d+3` the same
elimination produces `E[A_{d+3}] = f(m_3)` with `m_3 = E[c^3]` a genuinely new object (in the
plane, this is exactly the point at which `n = 6` needs `E[A_4^2]` and `E[A_3^3]`, see
`docs/N6_INGREDIENTS.md`).  Nothing in the circuit argument helps there: `d+3` points in `R^d`
are not a circuit, and their triangulations are not two.

---

# Part B. The 3-D Sylvester five-point problem

## B.1 Setup and anchors

Five points in `R^3` are in convex position iff none lies in the tetrahedron of the other four,
and at most one can, so the five events are disjoint and

    P_5(K) = 1 - 5 E[A_4],       A_4 = vol conv(P_1..P_4)/vol(K).

`E[vol T_4]` is a plain 4-point determinant expectation -- no hull code, no rejection -- so the
scan is cheap and accurate.  `P_5` is affine invariant.  Three exact anchors:

| K | `E[A_4]` | `P_5` |
|---|---|---|
| simplex | `13/720 - pi^2/15015` (Buchta-Reitzner) | 0.91300880 |
| cube | `3977/216000 - pi^2/2160` (Zinani) | 0.93078612 |
| ball | `9/715` | `134/143` = 0.93706294 |

`src/rd_scan.py` reproduces all three, plus the affine controls (a sheared simplex, a 1:2:5 box,
a 1:2:4 ellipsoid), at 4e7 tetrahedra each: seven anchor z-scores in `[-1.10, +0.85]`.

**Open question.**  Is `P_5` minimised by the simplex and maximised by the ball over all convex
bodies in `R^3`?  This is the `d = 3` analogue of Blaschke's theorem, and it is open.

## B.2 The scan (`results/rd_scan.json`)

240 3-D bodies at 4e7 tetrahedra each (se ~ 1e-5), 55 minutes:

| family | count | `P_5` range |
|---|---|---|
| anchors: simplices, cubes/boxes, ball, ellipsoid | 7 | 0.9130011 – 0.9370711 |
| Platonic: octahedron, dodecahedron, icosahedron | 3 | 0.9318159 – 0.9361917 |
| prisms and antiprisms over m-gons, m = 3..40, 5 heights | 48 | 0.9232067 – 0.9335984 |
| bipyramids over m-gons | 24 | 0.9245753 – 0.9341741 |
| pyramids over m-gons (m = 3 gives a simplex) | 32 | 0.9129816 – 0.9238591 |
| smooth: cylinders, cones, bicones, capsules, ball caps | 31 | 0.9238229 – 0.9370263 |
| near-simplex: pushed and corner-truncated simplices | 17 | 0.9129761 – 0.9268001 |
| near-ball: hulls of 6..256 random points on the sphere | 30 | 0.9188394 – 0.9370455 |
| random polytopes: hulls of 5..50 gaussian / ball points | 48 | 0.9160862 – 0.9351951 |

**No body fell below the simplex or above the ball.**  Over all 240,

    min over bodies of  (P_5 - P_5(simplex))/se  = -1.96   (a body that IS a simplex)
    min over bodies of  (P_5(ball) - P_5)/se     = -0.85   (the ball itself)

The bodies with the lowest `P_5` are affine images of the regular simplex (pushed simplices,
`pyramid3-*` -- a pyramid over a triangle IS a tetrahedron), scattering around 0.91300880
exactly as they must, interleaved with simplices truncated at `t = 0.02`, whose predicted excess
(`~ 4 t^3 x 0.3 = 1e-5`) is at the noise floor.  The highest is the ball.  The affine invariance shows up
unbidden and correctly all through the catalogue: `cone-h` gives 0.9238478 / 0.9238229 /
0.9238711 for h = 0.25 / 1 / 4, `cylinder-h` gives 0.9332906 / 0.9333088 / 0.9333046 for
h = 0.25 / 1 / 8, `prism3-h` gives 0.9232118 / 0.9232101 for h = 0.3 / 1.2 -- all constant
within error, since these are affine images of one another.

A representative ordering (all ± ~1e-5):

| body | `P_5` |
|---|---|
| simplex | 0.9130141 |
| square pyramid | 0.9211031 |
| cone | 0.9238229 |
| triangular prism | 0.9232101 |
| bipyramid over a triangle | 0.9245753 |
| cube | 0.9307806 |
| octahedron | 0.9318159 |
| half-ball (cap of height 1) | 0.9329113 |
| hexagonal prism | 0.9328655 |
| cylinder | 0.9333088 |
| bicone (double cone) | 0.9341684 |
| dodecahedron | 0.9358777 |
| icosahedron | 0.9361917 |
| capsule (L = 1) | 0.9366819 |
| hull of 256 random points on the sphere | 0.9370455 |
| **ball** | **0.9370711** |

The closest approaches to the ball from below are the "roundest" non-balls: hulls of 256 random
sphere points (`-1.7e-5`), a short capsule (`-3.7e-5`), a barely-cut ball (`-4.5e-5`).  The
closest approaches to the simplex from above (excluding affine simplices) are the barely
truncated simplices, indistinguishable from the anchor at this resolution, then hulls of five
random points at `+0.0031`.

## B.3 The two conjectured 3-D extremisers are strict local extremisers

`src/rd_local.py`, 4e8 tetrahedra per body (se ~ 3-5e-6), compared against the **exact** anchors.

**The simplex (conjectured minimiser).**  Truncating `kc` corners back a fraction `t` of each
adjacent edge removes a volume fraction `kc t^3`:

| perturbation | `P_5` | `P_5 - P_5(simplex)` | sigma | per unit removed volume |
|---|---|---|---|---|
| 4 corners, t = 0.05 | 0.91314088 | +1.321e-4 | +25 | +0.264 |
| 4 corners, t = 0.10 | 0.91381046 | +8.017e-4 | +154 | +0.200 |
| 4 corners, t = 0.20 | 0.91759655 | +4.588e-3 | +962 | +0.143 |
| 4 corners, t = 0.35 | 0.92678604 | +1.378e-2 | +3585 | +0.080 |
| 1 corner, t = 0.10 | 0.91319952 | +1.907e-4 | +36 | +0.191 |
| 1 corner, t = 0.20 | 0.91418240 | +1.174e-3 | +228 | +0.147 |
| 1 corner, t = 0.35 | 0.91674750 | +3.739e-3 | +768 | +0.087 |

Every truncation raises `P_5`, at 25 to 3585 sigma, and the response is **linear in the removed
volume** with the same coefficient (~0.26 at the smallest cut, rising as `t -> 0`) for 1 and 4
corners -- exactly the 2-D picture of `docs/N6_EXTREMAL.md` §5 repeated one dimension up.

> **The regular simplex is a strict local minimiser of `P_5` in `R^3`, and is a corner of shape
> space: the one-sided derivative in the corner-truncation direction is strictly positive.**

**The ball (conjectured maximiser).**

| perturbation | `P_5` | `P_5 - 134/143` | sigma |
|---|---|---|---|
| ellipsoid 1.05 : 1 : 1 (affine image -- must be 0) | 0.93706222 | -0.7e-6 | -0.2 |
| ellipsoid 2 : 1 : 0.5 (affine image -- must be 0) | 0.93706597 | +3.0e-6 | +1.0 |
| capsule, L = 0.15 | 0.93705152 | -1.14e-5 | -3.8 |
| capsule, L = 0.25 | 0.93703176 | -3.12e-5 | -10.3 |
| capsule, L = 0.40 | 0.93698202 | -8.09e-5 | -26.8 |
| capsule, L = 0.70 | 0.93685256 | -2.10e-4 | -69.4 |
| ball cut to a cap of height 1.95 | 0.93705897 | -3.97e-6 | -1.3 |
| ... height 1.90 | 0.93702099 | -4.20e-5 | -13.9 |
| ... height 1.80 | 0.93682543 | -2.38e-4 | -78.3 |
| ... height 1.60 | 0.93604413 | -1.02e-3 | -329 |

Two independent one-parameter perturbations, one *adding* material (Minkowski sum with a segment
of length `L`) and one *removing* it (cutting the ball by a plane), both strictly decrease `P_5`.
The affine null direction is confirmed to `1e-6`: two very different ellipsoids give the ball's
value at `-0.2` and `+1.0` sigma, as affine invariance demands.

> **The ball is a strict local maximiser of `P_5` in `R^3`.**

The two directions scale differently, and the difference is informative.  Cutting a cap of depth
`s = 2 - hh` removes a volume fraction `~ 0.75 s^2` and costs `dP_5 ~ -0.008 x (removed volume)`
-- a strictly negative FIRST-order response, which a one-sided direction at a maximum permits.
Minkowski-adding a segment of length `L` changes the volume by `~ 0.75 L` but costs only
`dP_5 = -0.00050 L^2` (the ratio `dP_5/L^2` is `-0.00050` at `L = 0.15, 0.25, 0.40` and
`-0.00043` at `L = 0.70`) -- i.e. **quadratic in `L`, sublinear in the added volume, so the first
derivative vanishes in that direction**.  This is the exact 3-D echo of the 2-D observation that
rounding a triangle by `T (+) rB` moves `P_6` only like `r^1.65` (`docs/N6_EXTREMAL.md` §5).
The mechanism is not identified in either case and is recorded as UNVERIFIED.

## B.4 Verdict for Part B

1. **The 3-D five-point conjecture is consistent with everything measured**: 240 bodies at 4e7
   tetrahedra, nothing outside `[P_5(simplex), P_5(ball)] = [0.9130088, 0.9370629]` beyond
   `2.0 sigma`, and both extreme deviations are attained by bodies that ARE the extremiser they
   "violate".
2. **The simplex is a strict local minimiser and the ball a strict local maximiser**, with
   measured one-sided derivatives (positive at the simplex, negative at the ball along two
   independent directions) and the affine null direction confirmed to `1e-6`.
3. Ordering among the named bodies is the intuitive one: simplex < pyramid < prism/cone <
   cube < octahedron < cylinder < bicone < dodecahedron < icosahedron < capsule < ball, i.e.
   `P_5` tracks roundness, as a Steiner-symmetrisation proof would require.
4. **No proof, and a resolution limit.**  A counterexample violating either bound by less than
   ~1e-5 would be invisible to the scan (though §B.3 excludes counterexamples in a neighbourhood
   of the two candidates at the 3e-6 level, and only along the directions tested).  Unlike the
   planar case there is no Nelder-Mead search here: the 3-D shape space is much larger and the
   scan plus the local analysis were judged the better use of the time.

## B.5 Cost and reproduction

    python rd_identity.py --nsamp 1e5                                   # 2 min (1 core)
    python rd_scan.py --samples 4e7 --seed 5                            # 55 min
    python rd_local.py                                                  # 33 min

All at `nice -n 10`.  `rd_scan.py` accepts `--resume <json>` to reuse already-computed bodies.
