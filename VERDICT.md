# VERDICT — convex-position probability, overnight server session 1

*2026-08-19, `~/math/convex-position-probability/`. Written against the brief in `TASK.md`.
Chronological detail is in `PROGRESS.md`; proofs in `THEOREMS.md`; sources and novelty in
`LITERATURE.md`; the n=5 table and scan in `docs/N5_LANDSCAPE.md`.*

## The short version

Three things were asked for. Two of them turned out to be **already solved in the
literature**, and finding that out is the session's most valuable result:

1. `P_disk(5) = 1 - 305/(48 pi^2)`, flagged "NOVELTY UNKNOWN" by the previous session, is
   **in Marckert 2017**, Table (7). Not new.
2. The `n = 5` extremal conjecture, which `TASK.md` set as the main overnight mission, is
   **a theorem of Marckert & Rahmani (2021)**. Only `n >= 6` is open.
3. `P_6` for the disk, `TASK.md`'s stretch goal, is **also in Marckert 2017** — and we were
   able to *extract* from his table the two moments the brief called "the hard part", and
   verify the whole thing closes exactly.

What genuinely came out of the night is (a) a new, fully verified exact tool (route P,
E1), (b) a clean structural theorem in the plane and in `R^3` that we could not find stated
anywhere (E2, E6), (c) an explicit closed form for `P_5` of every regular `m`-gon (E4), and
(d) a large amount of independent verification of the published record — including
reproducing Alikoski's 1939 formula to 50 digits and Buchta's entire expected-area table as
exact rationals, from scratch.

**Nothing here is claimed as new.** Identity (I) and its 3-D analogue are each one step from
three classical facts; identity (II) may be inside a paper we could not obtain. The novelty
column in `LITERATURE.md` says "not found", never "new".

## What is established

### E1. Route P — exact Sylvester quantities for any convex polygon *(new tool, verified)*

`src/polygon_exact.py`. Decompose the space of lines meeting a convex polygon into cells
indexed by the pair of edges a line crosses; parametrise a line in a cell by its two
boundary points. Santaló's `dp dtheta = sin(phi_i) sin(phi_j) ds_i ds_j / L` then makes the
`L^3` in `\int L^3 [c^k + (1-c)^k] dG` cancel exactly, and the integrand becomes a
**polynomial of degree `<= 2k+2`** in the two boundary parameters. Consequently

    T_k = J_k/(3|K|^2),   E[N_n] = C(n,2) T_{n-2},   E[A_{n-1}] = 1 - E[N_n]/n

are elementary and exact for every convex polygon — exact rationals for rational vertices
(`fractions.Fraction`), 50-digit for regular `m`-gons, and `O(1/m^4)`-convergent for smooth
bodies via inscribed polygons plus Richardson.

Evidence: two identities that must hold by construction (`T_0 = 2`, `E[N_3] = 3`) hold
exactly; it reproduces **Alikoski 1939** to 50 digits for `m = 3..24, 30, 40, 60, 100`; it
reproduces **Buchta 1984's entire published table** as exact rationals (square `11/144,
11/72, 79/360, 199/720`; triangle `1/12, 1/6, 43/180, 3/10, 197/560`); and it returns the
disk's `35/(48 pi^2)` to `6e-17`.

### E2. Identity (I) and its family *(proved; not found in the literature; do not claim new)*

For **every absolutely continuous probability measure on the plane** (no convexity, no
uniformity), with `A_k` the `mu`-content of the hull of `k` i.i.d. points:

    E[A_{n-1}] = 1 - (n-1) E[c^{n-2}],  c = mu-mass on one side of the line through two points,

and because `c` and `1-c` are identically distributed (swap the two points), the odd moments
of `c - 1/2` vanish, so a new free parameter enters only at every *other* step. Hence

    E[A_4] = 2 E[A_3]                              (I)
    E[A_6] = 3 E[A_5] - 5 E[A_3]                   (I')
    E[A_8] = 4 E[A_7] - 14 E[A_5] + 28 E[A_3]      (I'')  ... and so on.

Proof in `THEOREMS.md` §3 — three classical ingredients (Rényi–Sulanke, Efron, and a
one-line symmetry) and nothing else.

Evidence: exact `Fraction` arithmetic on six polygons (triangle, square, a sheared
parallelogram, the affinely-regular hexagon, a random rational pentagon, a random rational
heptagon) gives residual **exactly 0** for (I), (I') and (I''); sympy gives 0 identically on
the true regular `m`-gons `m = 3,4,5,6`; route P gives `< 1e-48` for `m = 3..60`; and Monte
Carlo confirms (I) for a standard Gaussian, a skewed anisotropic Gaussian mixture and a
heavy-tailed `t_3` law (`|z| < 1.7`).

**Novelty: NOT FOUND stated, but do not claim as new.** It is one step from three classical
facts, and it is *visible* in Buchta's published tables. The closest thing in the literature
is Finch ([arXiv:1601.04937](https://arxiv.org/pdf/1601.04937), p. 14) noting the Gaussian
instance as a curiosity — and that remark is about Lebesgue area, a logically different
statement from (I), which is about `mu`-content. (Both are true; both checked here.)

### E3. Identity (II) — `P_5` from two moments of the triangle area *(proved; status unresolved)*

    P_5 = 1 - 5 E[A_4] + 10 E[A_3^2]        (Efron-Buchta, classical)
        = 1 - 10 E[A_3] + 10 E[A_3^2]       (substitute (I))
        = 1 - 10 E[A_3 (1 - A_3)],

and for the uniform law on a convex body, `E[A_3^2] = (3/2) det Sigma_K / |K|^2`, so

    P_K(5) = (5/2) P_K(4) - 3/2 + 15 det(Sigma_K)/|K|^2 .

`P_5` depends on the body only through the first two moments of the random triangle's area:
Sylvester's functional and an affine-invariant covariance. **Status: UNRESOLVED.**
Marckert & Rahmani (2021) contain "a new formula for `Q_H^n` of independent interest" which
we could not read (Wiley HTTP 402; HAL behind an anti-bot wall). Until that paper is read,
presume (II) is known or nearly so.

### E4. Exact `P_5` for regular `m`-gons *(explicit closed form; not found tabulated)*

With `w = 2 pi/m`, from Alikoski + (II):

> **`P_5(regular m-gon) = 1 - 5 (15 cos^2 w + 92 cos w + 76) / (36 m^2 sin^2 w)`**

`m = 3` gives Valtr's `11/36`, `m = 4` gives Valtr's `49/144`, `m -> infinity` gives
`1 - 305/(48 pi^2)`. Sample values: `m=5: 7/12 - 47*sqrt5/450`, `m=6: 1373/3888`,
`m=8: 1469/2304 - 115*sqrt2/576`, `m=10: 461/720 - 229*sqrt5/1800`,
`m=12: 3439/5184 - 115*sqrt3/648`. Full table in `docs/N5_LANDSCAPE.md`.

### E5. `n = 6` for the disk, fully assembled and closed *(mostly from the literature)*

    P_6 = 1 - 6 E[A_5] + 15 E[A_4^2] - 20 E[A_3^3]

* `E[A_5](disk) = 7(2400 pi^2 - 3289)/(6912 pi^4)` — derived here (sympy, from
  `E[N_6] = 30 E[c^4]`), and confirmed by route P + Richardson to `2.9e-17`.
* `E[A_3^3](disk) = 1001/(6400 pi^4)` and `E[A_4^2](disk) = (2400 pi^2 + 31031)/(19200 pi^4)`
  — **extracted** from Marckert's `P^D_{6,3}` and `P^D_{6,4}` via `P(N_6=3) = 20 E[A_3^3]`
  and `P(N_6=4) = 15(E[A_4^2] - 4 E[A_3^3])`.

Cross-checks, all exact and all zero: Marckert's `P^D_{5,m}` and `P^D_{6,m}` each sum to 1;
our identities reproduce his `P^D_{5,3} = 15/(16 pi^2)` and `P^D_{5,4} = 65/(12 pi^2)`;
`E[N_6]` computed from his `P^D_{6,m}` equals our independent sympy value **exactly**; and
the three moments assemble to his `P_6 = 1 - (146400 pi^2 - 473473)/(11520 pi^4)` exactly.

### E6. A three-dimensional analogue *(proved; not found in the literature)*

The planar proof of (I) uses `#vertices = #edges`, which fails in `R^3`. But Euler's relation
plus `3F = 2E` gives `V = 2 + F/2` for a **simplicial** 3-polytope, and generic points give
one. With the `R^3` Rényi–Sulanke facet count and the same `c =_d 1-c` symmetry, Efron then
yields `E[A_4] = 3/5 - 2 m_2` and `E[A_5] = 3/2 - 5 m_2`, hence

> **`E[A_5] = (5/2) E[A_4]` in `R^3`, for every absolutely continuous law**
> (and `E[A_7] = (7/2) E[A_6] - (35/4) E[A_4]`).

So the expected volume of the hull of *five* points comes free from the *four*-point value:
`E[A_5](ball) = 9/286 = 0.0314685314...` and
`E[A_5](cube) = 3977/86400 - pi^2/864 = 0.0346069394...`.

Tested at 3e7 samples per body with `N_5`, `N_6` computed from exact orientation predicates
only (a point is a non-vertex iff it lies in `conv` of four of the others — Carathéodory —
so no hull construction and no floating-point volume enters): residuals `z = -1.30` (ball),
`-1.49` (cube), `+1.21` (simplex); the cube's `E[A_4]` hits Zinani's exact value at `z = +0.14`
and the ball's `E[A_5]` hits `9/286` at `z = +0.03`. `results/dim3_identity.json`.

Same novelty caveat as (I): one step from classical ingredients, not found stated, **do not
claim as new**.

## Evidence table

| claim | route 1 | route 2 | route 3 | agreement |
|---|---|---|---|---|
| `E[A_4] = 2E[A_3]` | exact `Fraction`, 6 polygons | sympy, regular `m`-gons `m=3..6` | 50-digit route P, `m=3..60`; MC on 3 non-uniform laws | exact 0 / `<1e-48` / `\|z\|<1.7` |
| `E[A_3]` regular `m`-gon | route P, 50 digits | Alikoski 1939 closed form | MC | `<1e-49` |
| `E[A_3..A_7]` square & triangle | route P exact rationals | Buchta 1984 published table | — | identical |
| `P_5(disk) = 1-305/(48pi^2)` | identity (II) + exact `E[A_3]`,`E[A_3^2]` | Marckert 2017 Table (7) | route P Richardson (`5.6e-16`); direct MC; route M | all agree |
| `P_disk(6,7,8)` | Marckert 2017 | our direct MC 2e8 | — | `\|z\| <= 1.55` |
| `E[A_5](disk)` | sympy `c^4` integral | route P Richardson | Marckert's `P^D_{6,m}` (exact) | `2.9e-17`, exact 0 |
| `P_5` on 9 bodies | identity (II) prediction | direct MC 2e8 | — | max `\|z\| = 2.47`, none above 3 |
| `P_5` Gaussian (non-uniform) | identity (II) from containment frequencies | direct MC 1.2e8 | classical `P_4 = (6/pi)arcsin(1/3)` | `z = -0.26`, `z = -0.60` |
| `n=5` extremal window | 561 exact bodies + optimisation | Marckert–Rahmani theorem | — | 0 violations |
| `E[A_5] = (5/2)E[A_4]` in `R^3` | Euler `V = 2 + F/2` proof | exact-predicate MC, 3 bodies | anchors 9/715 and Zinani | `\|z\| <= 1.5` |

## What failed / what was not achieved

* **No closed form for `E[A_5](regular m-gon)`.** The natural Alikoski-shaped ansatz
  `a_0 + P(cos w)/(m^2 sin^2 w) + Q(cos w)/(m^4 sin^4 w)` fitted at 60 dps on `m = 3..12`
  fails its held-out test at `m = 13..40` (residuals `2.7e-12` growing to `6.6e-10`, where a
  correct ansatz would give `~1e-45`). Reported as a failure.
* **No proof route for `n >= 6`.** Not attempted seriously once `n = 5` turned out to be a
  theorem; the honest note in `docs/N5_LANDSCAPE.md` §6 is that `F = E[A_3] - E[A_3^2]` is
  not obviously monotone under Steiner symmetrisation because both terms move the same way.
* **Marckert & Rahmani 2021 could not be obtained.** This is the one real gap: it decides the
  status of identity (II) and of E4.
* **Nelder–Mead under-converged for `<= 10`-gons** in the landscape optimisation (found
  `F = 0.0644406` where the regular decagon gives `0.0644200`). A limitation, not a result.
* **One bug** was found and fixed mid-session: a float64 covariance cancellation on
  near-degenerate polygons that made the optimiser report a spurious counterexample to the
  triangle bound (`PROGRESS.md` 00:15). **One methodological defect** was found and
  quantified: the A3 sweep reused one seed across bodies, correlating the z-scores
  (`PROGRESS.md` 01:10). **One estimator was rejected**: an early 3-D probe using scipy
  hull volumes suggested a 2.5σ excess for the ball that vanished under an exact-predicate
  estimator (`PROGRESS.md` 00:38).

### A3 in full: identities (I) and (II) on nine bodies

Reference `E[A_3]` exact (square, triangle), exact route P (pentagon, hexagon, octagon) or
route-P Richardson (disk, 3:1 ellipse, half-disk, stadium); `E[A_3^2]` exact throughout;
direct tester at 2e8 samples, moments at 1e8.

| body | reference for `E[A_3]` | z (I) | z (`P_4`) | z (II) | testers agree |
|---|---|---|---|---|---|
| square | exact rational | −0.06 | +1.51 | +0.40 | yes |
| triangle | exact rational | +1.68 | −0.32 | +2.09 | yes |
| disk | route P Richardson | +1.25 | −0.44 | +0.47 | yes |
| pentagon | route P exact polygon | +1.32 | +0.65 | −1.48 | yes |
| hexagon | route P exact polygon | +1.06 | +1.71 | −0.14 | yes |
| octagon | route P exact polygon | +1.16 | +0.03 | +0.96 | yes |
| 3:1 ellipse | route P Richardson | +1.25 | −0.44 | +0.47 | yes |
| half-disk | route P Richardson | −0.11 | +0.93 | +0.12 | yes |
| stadium (a=1) | route P Richardson | −1.04 | +1.33 | +2.47 | yes |

High-precision values for the two bodies with no literature value we found:
`half-disk: E[A_3] = 0.076512497523552, P_4 = 0.693950009905794, P_5 = 0.341069936631669`;
`stadium(a=1): E[A_3] = 0.074403975610728, P_4 = 0.702384097557086, P_5 = 0.352809265308243`.
Note the half-disk sits *below* the square in the four-point problem (`0.69395 < 25/36`) but
*above* it at `n = 5` (`0.34107 > 49/144`).

Two caveats, recorded rather than glossed:

* **The 3:1 ellipse is not an independent body.** Its sampler is the disk's scaled by
  `diag(3,1)` with the same seed, so it draws the affine image of the *same* point sets and
  reproduces the disk's z-scores digit for digit. That is a valid (and exactly passing)
  affine-invariance test of the pipeline, but it is one body, not two.
* **The 27 z-scores have mean +0.62**, which would be 3.2σ from zero *if they were
  independent*. They are not — every body used the same seeds. A dedicated check
  (`results/A3_bias_check.json`) settles it: square `P_4` against the exact `25/36` with
  **12 independent seeds** gives mean `z = -0.106 ± 0.289`, sd `1.065` (no bias, and the
  error bars are correctly sized); and re-running the whole body set with a fresh common
  seed moves the common offset from `+0.317` to `-0.193`. The stadium's `+2.47` becomes
  `-1.84` on a fresh seed. So the tilt carries no information about bias.

## Housekeeping notes for the operator

* `docs/PROBLEM_B.md` appeared in `docs/` at 23:45 (not written by this session). It briefs
  someone to attack the `n = 5` extremal conjecture **as an open problem**, which it is not.
  A clearly-marked correction banner has been prepended to it; the original text is intact
  below the banner. If another agent is working from that file, point them at the banner.
* `explainer/` is empty on the server and `README.md` still lists
  `explainer/in-convex-position.html`. That file exists only on the local machine; its DISK
  and regular-`m`-gon numbers are now exact and should be updated there before republishing.
* `PROGRESS.md` line 94 had a stray non-UTF-8 byte (from the previous session) that made
  `grep` treat the file as binary. Replaced with `-`.

## Recommended next mission

1. **Read Marckert & Rahmani, Mathematika 67 (2021) 860–884.** One paper decides whether
   E3 and E4 are new. Everything else in this campaign's "novelty" column is downstream of it.
2. If (II) turns out to be new or unstated, the natural write-up is short: identity (I) and
   its family (E2), (II) as a corollary, and the regular-`m`-gon table (E4) as the worked
   example — with route P as the computational appendix.
3. **Add a seed-variation rule to the harness.** Every sweep over configurations should
   vary the seed per configuration; test bias with many seeds on one configuration whose
   answer is exactly known. (See `PROGRESS.md` 01:10 — worth exporting to
   `../LESSONS_LEARNED.md`.)
4. **The open problem is `n >= 6`.** The n=6 analogue of (II) needs `E[A_4^2]` and
   `E[A_3^3]` for a general convex body. Route P does *not* reach them (they are not
   line-integral functionals). A concrete, self-contained target: find an exact route to
   `E[A_3^3]` (third absolute moment of the triangle determinant) for polygons — for the
   disk it is `1001/(6400 pi^4)`, which is suspiciously clean, and for the triangle and
   square it should be rational and may already be in Buchta's variance papers.
5. **3-D**: see E6 above — the analogue turned out to exist and be provable, and the
   obvious follow-ups are `R^d` for `d >= 4` (where `V = 2 + F/2` is replaced by the
   Dehn–Sommerville relations, so a collapse should survive only in odd dimensions) and a
   3-D analogue of identity (II).


---

## Addendum 2026-08-19 06:30 (local session, after the overnight run)

* **The n=5 proof route the night did not find now exists** (external model session prompted
  with `docs/PROBLEM_B.md`, relayed by the operator): a fibrewise section lemma showing that under
  Steiner symmetrisation and shaking the second triangle-area moment changes by exactly the
  fraction `b_2/V < 1` of the change in the first, so `F = E[A_3(1-A_3)]` is monotone and
  Blaschke's method applies. Verified here: the exact gap formulas reproduce
  `F(half-disk) - F(disk) = 0.001511838` to 9 digits against route P. Write-up:
  `docs/N5_PROOF.md`; `THEOREMS.md` section 8. The theorem stays attributed to
  Marckert-Rahmani 2021; the mechanism is recorded as an independent short proof.
* Elementary re-proofs of identities (I) and (II) (no Renyi-Sulanke / Efron / Buchta) are in
  the same place; (I) is now a one-line deterministic identity valid for every i.i.d. law.
* Item 1 of "Recommended next mission" (read Marckert-Rahmani 2021) is unchanged and still the
  gate on every novelty claim.

* **06:40 update: Marckert's paper obtained** as arXiv:1511.03658 (v1 2015; the Semantic Scholar
  DOI record points to it -- HAL/Wiley still block). Read in full: its "new formula" is a comb
  recursion (Props. 16-17) and its n = 5 argument is a polynomial comparison in a symmetry defect;
  identity (II), `E[A_4] = 2E[A_3]`, the covariance form and the regular-polygon values do NOT
  appear. Journal version (2021, with Rahmani) still unread. LITERATURE.md section 8.

* **n = 6 ingredient table (2026-08-19 morning):** E[A_3^3] and E[A_5] are now EXACT for
  triangle/square/disk via a Blaschke-Petkantschin width-function integral (E[A_3^k], any k,
  any convex body -- an extension of route P to area moments). E[A_3^3] = 31/9000, 137/72000,
  1001/(6400 pi^4). With E[A_5] = 43/180, 79/360 (Buchta) the full P_6 decomposition closes
  against Valtr (91/900, 49/400) and Marckert (disk), with E[A_4^2] MC-confirmed. See
  docs/N6_INGREDIENTS.md. Remaining open piece: an independent exact E[A_4^2 & convex] for a
  general body (convex-quadrilateral area second moment) -> first-principles P_6 for polygons.

---

## Addendum 2026-08-19 (server session, `TASK_N6.md`): E[A_4^2 & convex] for polygons — P_6 from first principles

**The last open piece of n = 6 is closed.** `src/n6_twochord_polygon.py` evaluates the two-chord
Blaschke-Petkantschin integral for `E[A_4^2 & convex]` on an arbitrary convex polygon. Method:
change variables from the two line offsets to the crossing point `X`
(`dp1 dp2 = |sin(phi1-phi2)| dX`); inside a cell of the arrangement cut by the lines through the
vertices parallel to `t(phi1)` and to `t(phi2)`, the four ray lengths are affine in `X`, so
`G1 G2` is a degree-10 polynomial and the inner 2-D integral is EXACT (degree-11 Duffy-Gauss on a
fan triangulation of the clipped cell). Only the 2-D angular integral is numerical: panelled at
every critical angle plus `phi2 = phi1`, graded, in 80-bit arithmetic.

* **Established.** `E[A_4^2 & convex](triangle) = 119/4500` to rel. 2e-20 and
  `(square) = 1307/54000` to rel. 1.3e-17, computed **without Valtr**; identical from unsymmetric
  vertex lists, from sheared affine images, and under raising the cubature order (which confirms
  the inner integral really is exact). So `P_6(triangle) = 91/900` and `P_6(square) = 49/400` are
  now **outputs** of this campaign's pipeline, and `P_6(K)` is assembled from first principles for
  any convex polygon.
* **New exact values** (regular m-gons, m = 5, 6, 8; ingredients and full table in
  `docs/N6_LANDSCAPE.md`, `results/n6_mgon_P6_final.json`):
  `P_6(pentagon) = 8941/22500 - 1349 sqrt5/11250`, `P_6(hexagon) = 461299/3499200`,
  `P_6(octagon) = 30103/61440 - 116141 sqrt2/460800`. 19-20 digit values for m = 7, 9, 10, 11, 12.
  Identification: smooth-denominator scan on the 80-bit values, calibrated by recovering
  119/4500 and 1307/54000. Cross-checked by direct convex-position Monte Carlo (6e8 samples).
* **Extremal conjecture (still open for n >= 6).** `P_6` is strictly increasing in `m` across
  the whole table, from 91/900 (triangle) to the disk value 0.1343093863571 — the ordering the
  conjecture predicts, now on exact values at m = 3, 4, 5, 6, 8 rather than Monte Carlo.
* **Negative result 1: there is no `P_5`-style closed form in `m` for n = 6.** Fitting
  `A(c)/S + B(c)/S^2` with `A,B in Q[c]`, `S = m^2 sin^2 w`, to 11 points good to 1e-18 gives a
  best 9-parameter residual of 8e-8, decaying smoothly with parameter count — interpolation, not
  discovery. A wider search over `Sum a_j c^j/(m^a sin^b w)` families found nothing for
  `E[A_3^3]` either, consistent with the earlier session's independent failure on
  `E[A_5](m-gon)`.
* **Negative result 2: the leading large-m correction is `1/m^4`, and its constant is not
  identified.** Exactly `P_5(disk) - P_5(m) = (7 pi^2/18)/m^4 + O(1/m^6)` (the `1/m^2` term
  cancels identically). For n = 6 the table was pushed to m = 20 and gives
  `a_4 = 3.2738126823 +- 1.3e-9`. A clean-looking PSLQ hit found at 6 digits,
  `a_4 = (4 pi^2 + 218 + 76/pi^2)/81` (one denominator, three small numerators, match 8e-10),
  was **refuted** by the sharpened value (3.6e-9 off, 2.9 spreads); the sharpened value yields a
  different "clean" relation at every basis/height. A textbook illustration of why the campaign
  requires PSLQ hits to come from many digits.
* **Bugs worth remembering.** Three precision bugs, none visible to grid refinement: a `float()`
  cast on the cell's affine coefficients inside the 80-bit path (biased every `F` by 8e-16
  relative), a decimal round-trip in the angle reconstruction, and float64 integration limits.
  Diagnosis required an independent exact-mpmath reference for the innermost quantity
  `F(phi1,phi2)`. Also: `((u+v)^5-u^5-v^5)/10` must be evaluated in its factored form
  `u v (u+v)(u^2+uv+v^2)/2` — the difference loses ~10 digits in sliver cells.

### Recommended next mission (updated)

1. Unchanged gate: read the journal version of Marckert-Rahmani 2021 before any novelty claim.
2. Push the m-gon table to m ~ 24 (the bottleneck is the mpmath width route for `E[A_3^3]`,
   `E[A_5]`, not the two-chord integral) and identify `a_4` for n = 6 exactly.
3. Get >= 25 digits on `E[A_4^2 & convex]` (double-double or a partly symbolic angular
   integration) to finish m = 10 and 12 exactly, and to reach m = 7, 9, 11 (cubic/quintic fields).
4. Run the n = 6 extremal scan over the 561-body family of `docs/N5_LANDSCAPE.md` — the machinery
   now exists for arbitrary polygons, so the n >= 6 conjecture can be tested at exact precision
   on non-regular bodies (half-disks, stadiums, one-vertex-pushed m-gons, random polygons).

* **Prior-art gate cleared (2026-08-19; docs/lit_triangle_area_moments.md).** E[A_3^3] for the
  triangle (31/9000) and square (137/72000) are PUBLISHED (Reed 1974; Beck 2024
  arXiv:2412.07952) -- NOT new, cite them. The genuine novelty candidates ("not found in
  accessible literature") are: the DISK higher triangle-area moments (MathWorld: distribution
  "apparently not known exactly"), E[A_4^2]/E[A_4^2&convex], the covariance-form P_5 identity, and
  the exact regular-m-gon P_5/P_6 (Morin 2024: "exact formulas are rare"). E[A_4]=2E[A_3] is
  known (Efron/Buchta). Residual unread source: Philip area12.pdf (404).

---

## Addendum (2026-08-19, server, `TASK_N6_EXTREMAL.md`): the n >= 6 extremal conjecture stress-tested, and the R^d identity

Full write-ups: [`docs/N6_EXTREMAL.md`](docs/N6_EXTREMAL.md), [`docs/RD_ANALOGUES.md`](docs/RD_ANALOGUES.md).

### Established (numerically, to the stated resolution)

* **The open conjecture `P_n(triangle) <= P_n(K) <= P_n(ellipse)` survives a serious attack at
  n = 6 and n = 7.**  568 area-normalised convex bodies (541 of them not affinely a triangle) at
  4e8 samples each: nothing outside the window beyond `1.3 sigma` at n = 6 or `2.2 sigma` at
  n = 7, and every extreme deviation is attained by a body that IS (an affine image of, or a
  polygonal approximation to) the extremiser it "violates".
* **The two conjectured planar extremisers are STRICT LOCAL extremisers, with measured shape
  derivatives.**  This is the part a scan cannot deliver.
  - Disk: the second variation along `h = 1 + eps cos(k theta)` is `-c_2 eps^2` with
    `c_2 = 0.3765(16), 0.5476(75), 0.6608(145), 0.7151(268), 0.676(63)` for `k = 3, 4, 5, 6, 8`
    (11 to 235 sigma), and `c_2(k=2) = +0.0005 +- 0.0004` -- the affine null direction, whose
    response is measured to be QUARTIC in `eps`, exactly as affine invariance requires.
  - Triangle: cutting corners raises `P_6` by `~0.30 x (removed area fraction)`, additively over
    corners; the one-sided derivative is strictly positive, so the triangle is a corner of shape
    space, not a smooth critical point.
* **Nelder-Mead over k-gons finds the triangle as the minimiser for every k = 4..8, and the
  REGULAR k-gon as the maximiser for every k = 4..8** (to ~1e-5).  Over support-function Fourier
  coefficients, the highest of ~6600 evaluations was `7.5e-5` below the disk.
* **Three apparent counterexamples were manufactured and killed by the protocol.**  A search
  started at the disk reported best-of-run gains of `+3.6e-5 ... +4.3e-5` over the disk (2.5-3
  sigma) -- a textbook selection maximum.  Re-measured at 1.2e10 samples with fresh seeds, all
  three sit AT or BELOW the exact disk value (`-1.08, -0.26, -1.21 sigma`).
* **The exact P_6 machinery now validated on NON-regular bodies.**  Seven irregular polygons
  (trapezoid, rhombus, truncated triangles, pushed pentagon, random quad/pentagon): the exact
  width + two-chord pipeline and 4e9-sample MC agree to `5e-6`, all `|z| <= 1.19`.  Seven new
  exact `P_6` values for non-regular bodies fall out.
* **`E[A_{d+2}] = ((d+2)/2) E[A_{d+1}]` in EVERY dimension `d >= 1`, POINTWISE.**  For `d+2`
  points in general position in `R^d`, `vol conv(all) = (1/2) sum_i vol conv(all minus i)`,
  because a circuit has exactly two triangulations (the two sides of its Radon partition).
  `d = 2` gives `E[A_4] = 2E[A_3]`, `d = 3` gives `E[A_5] = (5/2)E[A_4]` -- both previously
  derived here by Efron + Renyi-Sulanke + Euler.  Verified pointwise to `<= 7.1e-15` in
  `d = 2..7` and in expectation on 18 body/dimension pairs.
  **Not claimed as new** -- the circuit fact is standard, and `LITERATURE.md` already rules
  `E[A_4] = 2E[A_3]` known.  What is new is that one line replaces the whole chain and removes
  the dimension restriction.
* **The 3-D five-point conjecture (simplex min, ball max) is consistent**: 240 bodies at 4e7
  tetrahedra, nothing outside `[0.9130088, 0.9370629]` beyond `2.0 sigma`; and the simplex is a
  strict local minimiser (corner truncations raise `P_5` at 25-3585 sigma, linearly in the
  removed volume) while the ball is a strict local maximiser (both a Minkowski-added segment and
  a cut-off cap strictly lower `P_5`; ellipsoids give the ball's exact `134/143` to `1e-6`).

### Corrected / refuted

* **The R^3 note's conjecture that the "one relation collapses the vertex-count moment"
  phenomenon survives in ODD dimensions is WRONG in both directions.**  The Dehn-Sommerville
  route needs `f_{d-1} = alpha f_0 + beta` for simplicial `d`-polytopes, which holds iff
  `floor(d/2) = 1`, i.e. iff `d = 2` or `d = 3`; it fails for every `d >= 4`, even and odd alike
  (verified: in `R^4` `f_0 = 6` occurs with `f_3` in {8,9}; in `R^5` `f_0 = 7` with `f_4` in
  {10,12}; in `R^6` `f_0 = 8` with `f_5` in {12,15,16}).  The identity itself, by contrast, holds
  in every dimension.

### Recorded, unexplained

* At the triangle, Minkowski-rounding (`T (+) rB`) moves `P_6` only like `r^1.65` (fitted
  exponent; local slopes 1.5-1.9) instead of linearly.  A partial explanation: the INNER parallel
  bodies of a triangle are again triangles, so `r -> P_6(T (+) rB)` is constant on `r <= 0` and
  minimised there, forcing a vanishing one-sided derivative.  The exact exponent is UNVERIFIED.
* The 3-D echo: at the ball, `P_5(ball (+) segment_L) - P_5(ball) = -0.00050 L^2`, quadratic in
  `L` hence sublinear in the added volume, so that first derivative vanishes too.  Also
  UNVERIFIED as to mechanism.

### What would change the verdict

* At n = 6 a counterexample violating either bound by less than ~1e-5 is invisible to the scan
  (the local analysis of `docs/N6_EXTREMAL.md` §5 excludes them near the two candidates at a far
  finer level, but only along the directions tested).
* No search was run in 3-D shape space; only a 240-body scan plus a local analysis.  A
  Nelder-Mead / CMA search over 3-D polytope vertex coordinates is the obvious next step.
