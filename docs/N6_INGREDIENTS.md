# n = 6: the ingredient table, and an exact route to the triangle-area moments

*2026-08-19, local session. Code: `src/n6_bp_moments.py` (derivation + disk validation),
`src/n6_bp_polygon.py` (polygon width-function integrator), `src/n6_results.py` (assembly +
self-check), `src/mc_a3cube.py`, `src/mc_a4sq.py` (Monte-Carlo confirmations). Results JSON:
`results/n6_ingredients.json`.*

## The decomposition

    P_6 = 1 - 6 E[A_5] + 15 E[A_4^2] - 20 E[A_3^3],       (route_moments identity, n=6)

where A_k = (area of the convex hull of k uniform points)/|K|. Three moments enter. Two are
new relative to n <= 5:

* `E[A_3^3]` -- the third **absolute** moment of the triangle area. Unlike `E[A_3^2] = (3/2)
  det Sigma/|K|^2` (a polynomial moment), the odd power keeps the absolute value, so it is a
  genuinely new object.
* `E[A_4^2]` -- second moment of the 4-point hull area. It splits as
  `E[A_4^2] = E[A_4^2 & convex] + 4 E[A_3^3]`, the second term being the contribution of the
  configurations where one point is inside the triangle of the other three (then the hull is
  that triangle, and integrating the interior point gives a factor A_3, so `E[A_3^2 * A_3]`,
  times the 4 choices).

`E[A_5]` is not new: `E[A_5] = 1 - E[N_6]/6` with `E[N_6] = 15 E[(1-c)^4 + c^4]`, reachable by
the edge-count / area-fraction route.

## An exact route to E[A_3^k] for any convex body (Blaschke-Petkantschin width integral)

For three i.i.d. uniform points, `Delta = area = (1/2)|t1 - t2| d(P3, line(P1,P2))`. The planar
two-point Blaschke-Petkantschin identity `dP1 dP2 = |t1 - t2| dt1 dt2 dp dphi` (p = signed
distance of the line from the origin, phi = normal angle) separates the integral:

    E[A_3^k] = (1/V^{k+3}) (2^{1-k}/((k+2)(k+3)))
               int_0^pi int_p w(p,phi)^{k+3} [ int_t |t-p|^k w(t,phi) dt ] dp dphi,

where `w(.,phi)` is the width (slice-length) function of K in normal direction phi. This is
route P's philosophy (integrate over line space) applied to an area **moment** instead of an
edge count. For a convex polygon `w(.,phi)` is piecewise linear, so the inner t- and p-integrals
are polynomial on each piece and evaluated exactly by a fixed 8-node Gauss-Legendre rule; only
the phi-integral is adaptive, subdivided at every direction where two vertices project equally
(edges AND diagonals -- missing the diagonal directions was a real bug, caught by the square
failing at high precision).

Validation: reproduces the disk's `E[A_3^k]` for k=1,2,3 to > 30 digits
(`35/(48 pi^2)`, `3/(32 pi^2)`, `1001/(6400 pi^4)`), and the triangle/square `k=1,2` exact
rationals (`1/12, 1/72`; `11/144, 1/96`). New exact values (PSLQ from 40-digit integrals,
Monte-Carlo confirmed to rel. 1e-4):

> **E[A_3^3](triangle) = 31/9000,  E[A_3^3](square) = 137/72000,  E[A_3^3](disk) = 1001/(6400 pi^4).**

(NOVELTY, per docs/lit_triangle_area_moments.md: the triangle 31/9000 and square 137/72000
are KNOWN -- Reed 1974 / Beck 2024 arXiv:2412.07952 Tables 2.2/2.4 -- so cite them; the DISK value
1001/(6400 pi^4) is NOT FOUND in the literature -- MathWorld states the disk triangle-area
distribution is 'apparently not known exactly' -- and is a genuine novelty candidate.)

The same width machine gives `E[A_5]` independently (via `E[c^j]`, the area-fraction moments):
`E[A_5](triangle) = 43/180`, `E[A_5](square) = 79/360` -- matching Buchta 1984. And it confirms
`E[A_4] = 2 E[A_3]` to machine precision on every body.

## The verified n = 6 ingredient table

| body | E[A_5] | E[A_3^3] | E[A_4^2] | E[A_4^2 & convex] | P_6 | check |
|---|---|---|---|---|---|---|
| triangle | 43/180 | 31/9000 | 181/4500 | 119/4500 | **91/900** | Valtr |
| square | 79/360 | 137/72000 | 859/27000 | 1307/54000 | **49/400** | Valtr |
| disk | 7(2400 pi^2-3289)/(6912 pi^4) | 1001/(6400 pi^4) | (2400 pi^2+31031)/(19200 pi^4) | (2400 pi^2+19019)/(19200 pi^4) | **1 - (146400 pi^2-473473)/(11520 pi^4)** | Marckert |

Decimals: P_6 = 0.101111 (triangle), 0.122500 (square), 0.134309 (disk).

Provenance and independence:
* `E[A_5]`, `E[A_3^3]`: **exact**, from the B-P width route (independent of Valtr/Marckert).
* `E[A_4^2]`: for triangle/square, the value **forced** by the decomposition once `E[A_5]`,
  `E[A_3^3]` are known and `P_6` is set to Valtr's exact value; **independently confirmed by
  Monte Carlo** (`mc_a4sq.py`, 6e7 samples, rel. diff 1.6e-4 / 1.7e-5 / 1.3e-4). For the disk,
  `E[A_4^2]` is Marckert's, and the loop closes symbolically.
* The loop `P_6 = 1 - 6 E[A_5] + 15 E[A_4^2] - 20 E[A_3^3]` reproduces Valtr/Marckert exactly
  for all three bodies (`n6_results.py`: ALL LOOPS OK).

So for the triangle, square and disk the full n = 6 ingredient set is now pinned exactly and
cross-checked three ways (B-P, Valtr/Marckert, Monte Carlo).

## What is still open  (UPDATE: the method below now solves this; only the polygon coding remains)

**An independent exact route to `E[A_4^2 & convex]` (the convex-quadrilateral area second
moment) for a GENERAL convex body.** With it, `P_6(K)` would be assembled from first principles
for any K -- as `P_5` now is via `E[A_3(1-A_3)]` -- instead of leaning on Valtr/Marckert to
supply `E[A_4^2]`. A promising handle: for four points in convex position the hull area is
`(1/2)|det(diag_1, diag_2)|` (cross product of the two diagonals), which invites a two-line
Blaschke-Petkantschin. This is the concrete next target, and it is the last obstacle to a
closed-form `P_6` for regular polygons and other bodies (the disk already being done by
Marckert).

Once `E[A_4^2 & convex]` is in hand for polygons, `E[A_3^3]` already is (above), so
`P_6(regular m-gon)` would follow in closed form, extending the `P_5(regular m-gon)` result.


## The two-chord route to E[A_4^2 & convex] (the last piece) -- method validated

For four points in convex position the hull is a quadrilateral whose area is
`A_4 = (1/2)|det(diag_1, diag_2)|` (cross product of the two diagonals). Over the 3 ways to
split 4 labelled points into two pairs, exactly the crossing pair are the diagonals, so

    E[A_4^2 & convex] = (3/4) E[ D^2 1{P1P3 crosses P2P4} ],   D = det(P3-P1, P4-P2).

(Verified by Monte Carlo on triangle/square/disk: both sides match 119/4500, 1307/54000, and
(2400 pi^2+19019)/(19200 pi^4) -- `src/mc_diag_identity.py`.)

Now apply the two-point Blaschke-Petkantschin to EACH diagonal (line l1 through P1,P3 and line
l2 through P2,P4). With D = (s3-s1)(u4-u2) sin(phi1-phi2) and the two Jacobians |s1-s3|,
|u2-u4|, the point-integrals separate into a cubic straddle moment on each chord:

    E[A_4^2 & convex] = (3/4)(1/V^6) int_{l1} int_{l2} sin^2(phi1-phi2)
                          G(chord1; sigma1) G(chord2; sigma2) 1{X in K} dl1 dl2,

where X = l1 ^ l2, sigma_i is X's coordinate along line i, dl = dp dphi, and for a chord with
endpoints [a,b] along the line,
    G(a,b,sigma) = int int_{[a,b]^2} |s-s'|^3 1{s,s' straddle sigma} ds ds'
                 = ( (b-a)^5 - (b-sigma)^5 - (sigma-a)^5 ) / 10.
Note the 1/V^6 (not 1/V^4): D is a physical determinant, and E[(D/V)^2 ...] carries an extra
1/V^2 -- getting this wrong shows up as an exact factor of pi^2 on the disk.  Also G vanishes
at the chord endpoints, so the integrand is CONTINUOUS across the {X in K} boundary (it -> 0
there), which is why the quadrature converges fast.

Validation (disk, rotational symmetry collapses it to a 3-D integral over psi, p1, p2):
`src/n6_twochord_disk.py` reproduces (2400 pi^2 + 19019)/(19200 pi^4) to rel. 1.5e-12.

**Status: the method is validated; the general-polygon version is a 4-D line-space integral
(phi1, p1, phi2, p2) with polygon-line clipping for the chord endpoints and the crossing point.
Implementing it (server task) reproduces 119/4500 and 1307/54000 for triangle/square
independently of Valtr, and then yields closed-form P_6 for regular m-gons.**  With this, every
ingredient of P_6 is computed from first principles for an arbitrary convex body.

---

## RESOLVED (2026-08-19, server session, `TASK_N6.md`): the polygon two-chord integral

`src/n6_twochord_polygon.py` implements the two-chord integral for an **arbitrary convex
polygon**. Change variables from the two line offsets to the crossing point
(`dp1 dp2 = |sin(phi1-phi2)| dX`); inside a cell of the arrangement cut by the lines through the
vertices parallel to `t(phi1)` and to `t(phi2)`, the four ray lengths `u1,v1,u2,v2` are AFFINE in
`X`, so `G1 G2` is a degree-10 polynomial and the inner 2-D integral is done EXACTLY by a
degree-11 Duffy-Gauss rule on a fan triangulation of the clipped cell. Only the 2-D angular
integral is numerical (panelled at all critical angles, graded, 80-bit arithmetic).

Validation, independent of Valtr: `E[A_4^2 & convex](triangle) = 119/4500` to rel. 2e-20 and
`(square) = 1307/54000` to rel. 1.3e-17, with the same values recovered from unsymmetric vertex
lists and from sheared affine images, and bit-identical under raising the cubature order.
So **`P_6(triangle) = 91/900` and `P_6(square) = 49/400` are now OUTPUTS** of the pipeline.

New exact values for regular m-gons (full table, method, and the two negative results in
`docs/N6_LANDSCAPE.md`):

> `E[A_4^2 & convex]` = `1769/112500 + 577 sqrt5/168750` (pentagon), `403891/17496000` (hexagon),
> `35743/2764800 + 48793 sqrt2/6912000` (octagon), and hence
> **P_6(pentagon) = 8941/22500 - 1349 sqrt5/11250**, **P_6(hexagon) = 461299/3499200**,
> **P_6(octagon) = 30103/61440 - 116141 sqrt2/460800**.

Cross-checked by direct convex-position Monte Carlo (6e8 samples) and by `P_6` increasing
monotonically in `m` towards Marckert's disk value.

Three precision bugs had to be killed first, all of them invisible to grid refinement (the
lesson: a plateau that moves with neither the grid nor the arithmetic precision is a bug):
a `float()` cast on the cell's affine coefficients inside the 80-bit path (biased every `F` by
8e-16 relative), a decimal round-trip in the angle reconstruction, and float64 integration
limits. The innermost quantity `F(phi1,phi2)` had to be given an independent exact-mpmath
reference before the bug was findable.
