# The n = 6 landscape: P_6 for regular m-gons, from first principles

*Server session, 2026-08-19 (task `TASK_N6.md`). Code: `src/n6_twochord_polygon.py`
(the two-chord Blaschke-Petkantschin integral for `E[A_4^2 & convex]` on any convex polygon),
`src/n6_bp_polygon.py` (the width route for `E[A_3^3]` and `E[A_5]`),
`src/n6_mgon_final.py` (assembly), `src/pslq_n6.py` (exact identification),
`src/mc_mgon_n6.py` (independent Monte-Carlo check).
Results: `results/n6_mgon_P6_final.json`, `results/n6_tc_hiprec.json`,
`results/n6_mgon_mc_check.json`. Sibling of `docs/N5_LANDSCAPE.md`.*

## 0. What is new here

Before this session `P_6(K)` could be assembled only for the triangle, the square and the disk,
and only by **borrowing** `E[A_4^2]` from Valtr / Marckert (`docs/N6_INGREDIENTS.md`). The
missing piece was an independent route to the convex-quadrilateral area second moment
`E[A_4^2 & convex]`. It is now computed directly for an arbitrary convex polygon, so

> **every ingredient of `P_6` is computed from first principles for any convex polygon, and
> `P_6(triangle) = 91/900`, `P_6(square) = 49/400` are OUTPUTS of the pipeline, not inputs.**

New exact values (not found in the literature):

> **P_6(regular pentagon) = 8941/22500 − 1349√5/11250 = 0.129248382075802996**
> **P_6(regular hexagon)  = 461299/3499200 = 0.131829846822130773**
> **P_6(regular octagon)  = 30103/61440 − 116141√2/460800 = 0.133516325200578055**

together with their ingredients, and 19–20 digit values for m = 7, 9, 10, 11, 12.

## 1. The identity and the integral

`A_k` = (area of the convex hull of k uniform points)/|K|. With
`P_6 = 1 − 6E[A_5] + 15E[A_4^2] − 20E[A_3^3]` (the Buchta-type identity, `THEOREMS.md`) and
`E[A_4^2] = E[A_4^2 & convex] + 4E[A_3^3]`,

    P_6 = 1 − 6 E[A_5] + 15 E[A_4^2 & convex] + 40 E[A_3^3].

For four points in convex position the hull is a quadrilateral of area
`(1/2)|det(diag_1, diag_2)|`, and exactly one of the 3 pairings of 4 labelled points is the
crossing one, so `E[A_4^2 & convex] = (3/4) E[D^2 1{P1P3 crosses P2P4}]`, `D = det(P3−P1, P4−P2)`.
Applying the two-point Blaschke-Petkantschin identity to **each** diagonal separates the four
point integrals into one cubic straddle moment per chord (`docs/N6_INGREDIENTS.md`):

    E[A_4^2 & convex] = (3/4) V^{-6} ∫∫ sin^2(phi1−phi2) G(chord1;sigma1) G(chord2;sigma2)
                                        1{X in K} dp1 dphi1 dp2 dphi2,
    G(a,b,sigma) = ((b−a)^5 − (b−sigma)^5 − (sigma−a)^5)/10 .

## 2. How this file evaluates it (and why the value is trustworthy)

**Change of variables to the crossing point.** `dp1 dp2 = |sin(phi1−phi2)| dX`, so with
`u_i, v_i` = the two distances from `X` to `∂K` along `±t(phi_i)`,

    I = ∫_0^pi ∫_0^pi sin^2(D) |sin(D)| F(phi1,phi2) dphi1 dphi2,      D = phi1 − phi2,
    F(phi1,phi2) = ∫_K G1(X) G2(X) dX,   G_i = u_i v_i (u_i+v_i)(u_i^2+u_i v_i+v_i^2)/2,
    E[A_4^2 & convex] = (3/4) I / V^6 .

(The factored `G` is the same polynomial as `((u+v)^5−u^5−v^5)/10` but **cancellation-free**;
the naive difference loses ~10 digits in sliver cells where `u/v ~ 1e10`.)

**The inner 2-D integral is EXACT.** Cut `K` by the lines through each vertex parallel to
`t(phi1)`, and again by those parallel to `t(phi2)`. Inside one cell all four rays leave `K`
through fixed edges, so `u1,v1,u2,v2` are **affine in X** and `G1·G2` is a polynomial of total
degree 10. A degree-11-exact 6×6 Duffy-Gauss rule on a fan triangulation of the (Sutherland-
Hodgman clipped) cell therefore integrates it exactly. Confirmed: raising the rule order to
ng = 8 or 10 changes nothing (bit-identical in float64; differences <= 1e-21, i.e. rounding, in
80-bit).

**Only the 2-D angular integral is numerical.** It is panelled at every critical angle (the
directions of all vertex differences and their perpendiculars — where the strip combinatorics
changes and `t(phi)` becomes parallel to an edge) and at `phi2 = phi1` (the `|sin D|^3` kink),
with graded sub-panels, then Gauss-Legendre on each panel. Every factor in the integrand is
`>= 0`, so the sum has no cancellation. For a body invariant under rotation by `2pi/m`, `phi1`
is restricted to `[0, pi/q)` with multiplicity `q = m` (m odd) / `m/2` (m even) — the order of
the induced shift group on `phi` in `R/pi Z`.

**Arithmetic.** 80-bit (`np.longdouble`) throughout, with Gauss nodes, `pi`, the integration
limits and all trigonometric values generated in mpmath and cast down. This is what buys the
last four digits: in double precision the same code saturates near 1e-15 relative (and near
1e-10 if the angular panels are not graded).

**Validation (independent of Valtr).**

| check | result |
|---|---|
| `E[A_4^2 & convex](triangle)` vs `119/4500` | rel. 2.1e-20 (best of 4 settings), 3.2e-19 (worst) |
| `E[A_4^2 & convex](square)` vs `1307/54000` | rel. 1.3e-17 |
| same values from the *unsymmetric* vertex lists (q = 1) and from sheared affine images | agree |
| rule order ng = 6 / 8 / 10 | identical to rounding (inner integral is exact) |
| `P_6(triangle)` assembled | 0.101111111111111111239 vs `91/900` (1.3e-19) |
| `P_6(square)` assembled | 0.122499999999999986519 vs `49/400` (1.3e-17) |
| `m → ∞` | 0.1341519 (m=12) → disk 0.1343093863571 (Marckert) |
| direct convex-position Monte Carlo, 6e8 samples each (`mc_mgon_n6.py`) | pentagon z = +0.00, hexagon z = +1.39, octagon z = +1.37 |

## 3. The exact table

`E[A_5]` and `E[A_3^3]` are from the width route (`n6_bp_polygon`), PSLQ'd from 20+-digit
integrals; `E[A_4^2 & convex]` is from the two-chord integral above, identified by a
smooth-denominator scan on the 19–20 digit value (calibrated: the scan returns `119/4500` and
`1307/54000` for m = 3, 4). `w = 2 pi/m`, `c = cos w`.

| m | E[A_5] | E[A_3^3] | E[A_4^2 & convex] | E[A_4^2] | **P_6** | P_6 decimal |
|---|---|---|---|---|---|---|
| 3 | 43/180 | 31/9000 | 119/4500 | 181/4500 | **91/900** | 0.101111111111111111 |
| 4 | 79/360 | 137/72000 | 1307/54000 | 859/27000 | **49/400** | 0.122500000000000000 |
| 5 | 221/1500 + 34√5/1125 | 32/28125 + 19√5/75000 | 1769/112500 + 577√5/168750 | 2281/112500 + 374√5/84375 | **8941/22500 − 1349√5/11250** | 0.129248382075802996 |
| 6 | 149347/699840 | 57709/34992000 | 403891/17496000 | 57701/1944000 | **461299/3499200** | 0.131829846822130773 |
| 7 | 0.212777561535668250 | 0.001628085431817641 | 0.022966895063629595 | 0.029479236790900159 | — | 0.132961474013140063 |
| 8 | 4531/36864 + 5837√2/92160 | 193/230400 + 53√2/96000 | 35743/2764800 + 48793√2/6912000 | 45007/2764800 + 64057√2/6912000 | **30103/61440 − 116141√2/460800** | 0.133516325200578055 |
| 9 | 0.212325614582656606 | 0.001613490463064818 | 0.022881818427430083 | 0.029335780279689355 | — | 0.133813207438104336 |
| 10 | 361891/3000000 + 5121√5/125000 | 92369/112500000 + 52973√5/150000000 | 0.022865299252066932 | 0.029308229802433946 | — | 0.133983396871405807 |
| 11 | 0.212184680329545883 | 0.001609093492942185 | 0.022855391229180319 | 0.029291765200949060 | — | 0.134086526178116899 |
| 12 | 253843/2239488 + 79843√3/1399680 | 6103/7464960 + 1183√3/2592000 | 0.022849154647336518 | 0.029281428301003499 | — | 0.134151930266256898 |
| 14 | 0.212114697948167428 | 0.001606949269692774 | 0.022842302788476966 | 0.029270099867248063 | — | 0.134224324925860894 |
| 16 | 0.212096995434246596 | 0.001606412678334944 | 0.022838997718677478 | 0.029264648432017253 | — | 0.134259500308080344 |
| 18 | 0.212087605128454468 | 0.001606129362967592 | 0.022837245939543922 | 0.029261763391414289 | — | 0.134278232841135689 |
| 20 | 0.212082251266260088 | 0.001605968344592022 | 0.022836247727072553 | 0.029260121105440641 | — | 0.134288942092208652 |
| ∞ (disk) | 7(2400π²−3289)/(6912π⁴) | 1001/(6400π⁴) | (2400π²+19019)/(19200π⁴) | (2400π²+31031)/(19200π⁴) | **1 − (146400π²−473473)/(11520π⁴)** | 0.134309386357109939 |

(Rows m = 14..20 were computed for the asymptotics in section 5(b) with a single (nphi, grade)
setting rather than three, so their error bar is ~1e-19 rather than ~1e-20.)

For m = 7, 9, 11 the field `Q(cos 2pi/m)` has degree 3, 3, 5 over `Q` and the exact value is out
of PSLQ reach at 19 digits (a plain PSLQ run returns relations with denominators like
340946 / 415187 / 2513 — noise, and reported as such). For m = 10 and 12 the field is quadratic
but the denominators are evidently larger than `3e6 x` (nothing found for smooth `D <= 1e9`);
these need >= 25 digits, which the 80-bit path cannot deliver.

## 4. The extremal conjecture at n = 6

`P_6` is **strictly increasing in m** over the whole table and converges to the disk value:

    0.1011111 (m=3) < 0.1225000 < 0.1292484 < 0.1318298 < 0.1329615 < 0.1335163
        < 0.1338132 < 0.1339834 < 0.1340865 < 0.1341519 (m=12) < 0.1343094 (disk).

This is exactly the ordering the (still **open**) extremal conjecture
`P_triangle(n) <= P_K(n) <= P_ellipse(n)` predicts for n = 6 — Marckert-Rahmani 2021 proved only
n = 5. The regular m-gon family is of course a one-parameter family and monotonicity along it is
weak evidence; but it now rests on **exact** values at m = 3, 4, 5, 6, 8 rather than Monte Carlo.

## 5. Two negative results, recorded honestly

**(a) There is no `P_5`-style closed form in `m`.** For n = 5 everything collapses to
`P_5 = 1 − 5(15c²+92c+76)/(36 m² sin²w)` — a single polynomial in `c` over `m² sin²w`. The
n = 6 analogue **fails**. Fitting `X = A(c)/S + B(c)/S²` with `A, B in Q[c]`, `S = m² sin²w`, to
the 11 data points (m = 3..12 and the disk as `c = 1, S = 4π²`), all of which are good to
~1e-18, the best 9-parameter residual is **8e-8** and the residual decays smoothly with the
number of parameters instead of dropping to the noise floor — i.e. the fits are interpolating,
not discovering a form. A wider automated search over families
`X = Σ_j a_j c^j /(m^a sin^b w)` (a, b <= 8, up to 3 groups, <= 8-10 parameters) found nothing for
`E[A_3^3]` either. This is consistent with the earlier session's independent failure to fit
`E[A_5](regular m-gon)` (`results/C_mgon_EA5_fit.json`, verdict "fails"). The individual values
*are* algebraic numbers in `Q(cos 2pi/m)` (see the table) — there just seems to be no uniform
rational-function-of-`c` formula.

**(b) The leading large-m correction is `1/m^4`, not `1/m^2`, and its constant resisted
identification.** For n = 5 the closed form gives exactly

    P_5(disk) − P_5(regular m-gon) = (7 pi^2 / 18) / m^4 + O(1/m^6),      7 pi^2/18 = 3.8381795,

the `1/m^2` term cancelling identically. The same cancellation happens at n = 6:
`(P_6(disk) − P_6(m)) m^4` runs 2.6891, 3.0232, 3.1631, 3.2135, 3.2363, 3.2484, 3.2554, 3.2599,
3.2629, 3.2650 (m = 3..12), 3.26772, 3.26933, 3.27037, 3.27108 (m = 14, 16, 18, 20). The table was
extended to m = 20 for this purpose (`results/n6_mgon_large_m.json`; the mpmath width route for
`E[A_3^3]`, `E[A_5]` costs 45 min at m = 20, the two-chord integral a few minutes). Richardson
fits in `1/m^2` with 5-10 terms then agree to 1.3e-9:

    a_4 = 3.2738126823 +- 1.3e-9      (a_4 / pi^2 = 0.3317065760, vs 7/18 = 0.3888889 for n = 5).

**And that is as far as it goes: `a_4` is NOT identified.** At 6 digits (m <= 12) PSLQ over
`{1, pi^2, pi^-2}` returned the strikingly clean `a_4 = (4 pi^2 + 218 + 76/pi^2)/81` — one
relation, one denominator, three small numerators, matching to 8e-10. Extending the table to
m = 20 moved `a_4` by 3.6e-9 and **refuted it** (the candidate now sits 2.9 spreads away), and the
sharpened value produces different "clean" hits at every basis and coefficient bound
(`11/36, 19/180, 271/180` at height 300; `4/7, -57/14, 247/14, -8` with a `pi^-4` term). All noise.
Nine digits is simply not enough for a 3-4 term relation; a proper identification needs the
asymptotic expansion done symbolically, or `a_4` to ~1e-15.

## 6. Cost

The two-chord integral at 19-digit accuracy costs 20 s (m = 3) to 12 min (m = 11) on 14 cores.
The bottleneck for the table as a whole is the mpmath width route for `E[A_3^3]` and `E[A_5]`.
