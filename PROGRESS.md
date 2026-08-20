# PROGRESS — convex-position probability campaign

Running log, `date '+%F %T'` timestamps (Europe/Prague local). Audience-facing
status lives in README.md / VERDICT.md; this file is the chronological record.

## 2026-08-18 22:52 — campaign opened
Read `Math/better-research-targets-than-3x3-eigenvalues.md` (Candidate 2 =
Sylvester/Valtr convex-position probabilities) and `Math/LESSONS_LEARNED.md`.
Folder created locally (`Math/convex-position-probability/`) and on the compute
server (`<server>:~/math/convex-position-probability/`, venv with
numpy 2.5.2 / scipy 1.18 / sympy 1.14 / mpmath 1.3 / numba 0.67).
Server has an unrelated live Claude session in `~/math/nonmonic-01` and a Lean
build in `~/math/nonmonic-cubic-lean` — do not touch either.

## 2026-08-18 22:53 — 8-angle literature sweep launched (workflow wf_cb73e498)
Angles: classical Sylvester; Valtr exact; disk/Marckert; Barany asymptotics +
extremal conjecture; higher dims; Gaussian/recent arXiv; computational/OEIS;
open-problem surveys. Then completeness critic, 6 adversarial verifiers on the
load-bearing formulas/values, and a synthesis into LITERATURE.md.

## 2026-08-18 22:55 — anchors + MC pipeline written, validated locally
`src/anchors.py`: Valtr parallelogram/triangle formulas (Fraction-exact),
Sylvester disk 1-35/(12 pi^2), 3D five-point cube (Zinani) and ball (9/715).
Sanity: n=4 gives 25/36, 2/3; n=5 square 49/144.
NOTE: the targets doc quotes the disk four-point value as "~0.7037"; the correct
decimal of 1-35/(12 pi^2) is 0.704479881... (typo in the doc, not in the math).
`src/convex_position.py`: two structurally independent testers (numpy
triangle/tetrahedron containment vs numba monotone-chain hull) that must agree
sample-for-sample; samplers for square/disk/triangle/regular k-gon/cube/ball/
simplex. `--validate` reproduces all 9 anchors: local 4e5-sample run all |z|<2,
testers agree exactly (pure-python fallback for tester B locally, no numba).

## 2026-08-18 22:57 — server: 2e7-sample anchor validation launched (nice -n 10)
`results/validation_anchors_2e7.{log,json}`.

## 2026-08-18 23:03 — server anchor validation PASS 9/9 at 2e7
All |z| < 3.1 (square_n5 z=-3.02 is the outlier; queued a 4e8 recheck with a
fresh seed and a 5e7 both-testers run). Testers agree exactly on every anchor.
`results/validation_anchors_2e7.json`.

## 2026-08-18 23:06 — first MC sweep done (results/sweep.log, *_2e8.json)
Disk, 2e8 samples each (numba tester B; A/B cross-check on 2e6 agrees exactly):
  n=5  0.35618818 +- 3.4e-5
  n=6  0.13432961 +- 2.4e-5
  n=7  0.0391123  +- 1.4e-5
  n=8  0.00904229 +- 6.7e-6
  n=9  0.001698025 +- 2.9e-6
  n=10 0.000262035 +- 1.1e-6
Square/triangle n=7..10 at 1e8 vs Valtr exact: all |z| <= 2.13 (8 anchors).
Ordering disk > square > triangle holds at every n (consistent with the
P_triangle <= P_K <= P_disk conjecture). These disk numbers are UNVERIFIED
against the literature until Marckert's recursion values are compared.

## 2026-08-18 23:10 — Route M (hull-area moments) written and validated
`src/route_moments.py`: P_n = sum_{j=0}^{n-3} (-1)^j C(n,j) E[A_{n-j}^j]
(Buchta-type identity, derived from P(N_n=m) = C(n,m) E[1_conv(m) A_m^{n-m}];
exact rational coefficients printed by `coefficients(n)`; e.g.
P_5 = 1 - 5E[A_4] + 10E[A_3^2]). Jarvis-march hull area (no code shared with
the direct testers). Server validation 4e6 samples: 7/7 planar anchors, |z|<1.5.
Disk n=5,6,7 route-M runs at 4e8 launched (results/disk_n*_routeM_4e8.json).

## 2026-08-18 23:20 — Route S: semi-analytic disk n=4,5 -> P_5(disk) EXACT
`src/disk_semi_analytic.py` (mpmath, 40 digits) + `src/disk_n5_exact.py` (sympy).
Ingredients: Renyi-Sulanke edge count E[N_n] = C(n,2) E[(1-c)^{n-2} + c^{n-2}]
over the chord through two random points (Blaschke-Petkantschin: ordered pair
gives int int |t1-t2| = L^3/3 -- NOTE first attempt used L^3/6 and was off by
exactly 2x; caught by the P_4 anchor, which then reproduced to 2e-41);
Efron E[A_{n-1}] = 1 - E[N_n]/n; Buchta P_5 = 1 - 5E[A_4] + 10E[A_3^2];
E[A_3^2] = 3/(32 pi^2) exactly (E[det^2] = 3/8, polynomial moment).
Because (1-c)^3 + c^3 = 1 - 3c + 3c^2, n=5 needs only the same integrals as n=4.
RESULTS (sympy exact, p = cos(phi) substitution makes all integrals elementary):
    E[N_4](disk) = 4 - 35/(12 pi^2)          -> P_4 = 1 - 35/(12 pi^2)  [anchor OK]
    E[N_5](disk) = 5 - 175/(24 pi^2)
    E[A_4](disk) = 35/(24 pi^2)  = 2 * E[A_3]   (curious: exactly twice)
    P_5(disk)    = 1 - 305/(48 pi^2) = 0.35618831227264541...
Cross-checks: direct MC 2e8: 0.35618818 +- 3.4e-5 (z = -0.4); PSLQ on the
40-digit numeric found [48,-48,305,0] unprompted (basis 1, pi^-2, pi^-4, pi^-6).
Route M 4e8 pending. NOVELTY UNKNOWN until LITERATURE.md lands (Marckert 2017
computes disk values; this may well be there) -- do not claim new.
Consequence for n=6: (1-c)^4 + c^4 does not collapse, and P_6 also needs
E[A_4^2], E[A_3^3] which the edge-count route does not give -> harder.

## 2026-08-18 23:16 — cross-checks in, session paused (operator restarting PC)
Route M disk n=5 (4e8/k): 0.35618671 +- 1.7e-5 vs exact 1-305/(48pi^2) = 0.35618831 (z=-0.09).
  Moments: E[A_4] = 0.14776048 (exact 35/(24pi^2) = 0.14776006), E[A_3^2] = 0.0094989 (exact 3/(32pi^2) = 0.0094986).
  => P_5(disk) now confirmed THREE independent ways (direct MC, route M, exact route S).
Route M disk n=6: 0.13432529 +- 3.9e-5 vs direct 0.13432961 +- 2.4e-5 (z~0.1). Agree.
Square n=5 recheck: 4e8 -> 0.3402545 +- 2.4e-5 (z=-0.98 vs 49/144); 5e7 both testers 0.3403399 (z=+0.93, agree).
  The earlier z=-3.02 was a fluctuation. Anchor set fully clean.
Still running on server (nohup, nice): route M disk n=7 (results/disk_n7_routeM_4e8.json).
Literature workflow (wf_cb73e498) was still in its 8-agent sweep stage when the
session paused; see RESUME.md for how to recover its transcripts.

## 2026-08-18 23:19 - last in-flight server job finished
Route M disk n=7 (4e8/k): 0.0390578 +- 5.9e-05 vs direct 0.0391123 +- 1.4e-05 (z=-0.90). Agree.
No server jobs from this campaign remain running.

## 2026-08-18 23:36 — session resumed after PC restart; handing off to a server tmux session
General identities derived (see TASK.md, to be proved/verified by the server session):
  (I)  E[A_4] = 2 E[A_3] for every planar convex body (edge counts + E[c]=1/2 + Efron);
  (II) P_5(K) = 1 - 10 E[A_3(1-A_3)] = (5/2)P_4 - 3/2 + 15 det(Sigma_K)/|K|^2
       (checks: triangle 11/36 w/ E[A_3^2]=1/72; square 49/144 w/ 1/96; disk 1-305/(48pi^2)).
Local literature agents (4, in background) write to docs/lit_*.md; server session does
its own novelty check meanwhile (TASK.md section A4).

## 2026-08-18 23:47 — [server session] A1 done: all three code paths reproduce the anchors
Overnight server session started 23:40 (tmux `convex-position`). Re-ran every path first:
- `convex_position.py --validate` at 2e7: 9/9 anchors, |z| <= 1.58, testers A and B agree
  exactly on every anchor (`results/A1_validate_direct_2e7.{log,json}`).
- `route_moments.py --validate` at 4e6: 7/7 planar anchors, |z| <= 1.47
  (`results/A1_validate_routeM_4e6.{log,json}`).
- `general_n5_identity.py`: E[A_4]/E[A_3] = 1.997-2.000 on square/triangle/disk/pentagon/hexagon
  (2e6 triples each, consistent with identity (I)); identity (II) reproduces 49/144, 11/36 and
  1 - 305/(48 pi^2) to all printed digits (`results/A1_general_n5_identity.log`).
Note: `docs/` is EMPTY on the server — the local session's `docs/lit_*.md` were never synced,
so this session does its own literature sweep (TASK A4).

## 2026-08-18 23:55 — [server] NEW TOOL "route P": exact Sylvester quantities for ANY convex polygon
`src/polygon_exact.py`. Idea: decompose the space of lines meeting a convex polygon into
cells indexed by the unordered pair {i,j} of edges the line crosses, and parametrise a line
in a cell by its two boundary points P(u) = A_i + u d_i, Q(v) = A_j + v d_j.  Santalo's
  dp dtheta = sin(phi_i) sin(phi_j) ds_i ds_j / L
becomes  dG = |d_i x (Q-P)| |d_j x (Q-P)| / L^3 du dv, so in
  J_k = int L^3 [c^k + (1-c)^k] dG        (c = cut-off area fraction)
the L^3 CANCELS EXACTLY and the integrand is a POLYNOMIAL of degree <= 2k+2 in (u,v):
  J_k = sum_{i<j} int_0^1 int_0^1 (-(d_i x (Q-A_i))(d_j x (A_j-P))) [c^k+(1-c)^k] du dv,
c(u,v) = shoelace([P, V_{i+1},...,V_j, Q])/|K|, degree <= 2.  Hence
  T_k = J_k/(3|K|^2),  E[N_n] = C(n,2) T_{n-2},  E[A_{n-1}] = 1 - E[N_n]/n
are ELEMENTARY and EXACT for every convex polygon.  No quadrature error: Gauss-Legendre
with enough nodes is exact for polynomials; a pure `fractions.Fraction` implementation
(`J_integrals_exact_rational`) gives exact rationals for rational vertices.
Built-in checks that must hold identically: T_0 = 2 (i.e. int L^3 dG = 3|K|^2) and
E[N_3] = 3 T_1 = 3.  Both hold exactly (Fractions) and to 1e-50 (mpmath).

## 2026-08-18 23:56 — [server] A3/A2 evidence: identity (I) E[A_4] = 2 E[A_3] verified EXACTLY
`src/exact_rational_check.py` (exact Fractions, no floating point at any step):
  body                       E[A_3]        E[A_4]-2E[A_3]   P_4          E[A_3^2]   P_5
  triangle                   1/12          0 (exact)        2/3          1/72       11/36
  unit square                11/144        0 (exact)        25/36        1/96       49/144
  sheared parallelogram      11/144        0 (exact)        25/36        1/96       49/144
  affinely-regular hexagon   289/3888      0 (exact)        683/972      25/2592    1373/3888
  random rational pentagon   68598907/915546564    0 (exact)  ...        ...        ...
  random rational heptagon   17956351/241813452    0 (exact)  ...        ...        ...
So (I) is not an artefact: it holds as an identity of exact rationals on arbitrary polygons.
Every literature anchor (Valtr 2/3, 25/36, 11/36, 49/144; hexagon 683/972) is reproduced
exactly. Bonus exact values from the same run (E[A_5], first ingredient of P_6):
  E[A_5](triangle) = 43/180,  E[A_5](square) = 79/360,  E[A_5](aff.-reg. hexagon) = 149347/699840.
Also exact via sympy on the true regular m-gons (irrational vertices): m=3,4,5,6 all give
E[A_4] - 2E[A_3] = 0 identically (`mgon_exact.py exact`, m=5 took 397 s).

## 2026-08-19 00:13 — [server] A4 LITERATURE SWEEP DONE — two "results" turn out to be KNOWN
`LITERATURE.md` written. Headline verdicts, both from primary sources read in full:
1. **P_5(disk) = 1 - 305/(48 pi^2) IS IN THE LITERATURE** — Marckert 2017 (arXiv:1402.3512),
   Table (7): `1 - P^D_5 = 305/(48 pi^2)`. The local session's "NOVELTY UNKNOWN" resolves to
   KNOWN. The same table gives exact P_6, P_7, P_8 for the disk; our MC agrees with all of
   them (|z| <= 1.55), which validates the pipeline and confirms the table:
     n=6 exact 0.1343093864 vs MC 0.13432961 +- 2.4e-5 (z=+0.84)
     n=7 exact 0.0390905623 vs MC 0.0391123  +- 1.4e-5 (z=+1.55)
     n=8 exact 0.0090342622 vs MC 0.00904229 +- 6.7e-6 (z=+1.20)
   => TASK.md section C ("P_6 for the disk") is a solved problem; do not spend the night on it.
2. **The n=5 extremal conjecture IS A THEOREM** — Marckert & Rahmani, *Around Sylvester's
   question in the plane*, Mathematika 67 (2021) 860-884, proved P_T(5) <= P_K(5) <= P_D(5).
   Verified by a verbatim quote in Marckert & Morin, arXiv:2411.08456 p.4 (read in full):
   "The case n = 5 was proved by Marckert and Rahmani [18], but n >= 6 is still a conjecture".
   => TASK.md section B3 is also a solved problem; our landscape becomes a VERIFICATION.
   Their paper contains "a new formula for Q_H^n of independent interest" which we could NOT
   read (Wiley 402, HAL anti-bot). That formula may well be identity (II). Recorded as the
   single most important gap; identity (II) is therefore logged as UNRESOLVED, presumed known.
3. Alikoski 1939's regular-m-gon mean triangle area (via MathWorld "Polygon Triangle Picking"),
   E[A_3] = (9cos^2 w + 52 cos w + 44)/(36 m^2 sin^2 w), w = 2pi/m, is INDEPENDENTLY REPRODUCED
   by route P: exactly for m=3,4,6, and to 50 digits for m = 3..24,30,40,60,100 (max residual
   3.4e-49). MathWorld's tabulated values (9+2sqrt5)/180 (m=5) and (97+52sqrt2)/2304 (m=8) match.
4. Route P also reproduces **Buchta's entire published expected-hull-area table** exactly, as
   rationals: square 11/144, 11/72, 79/360, 199/720; triangle 1/12, 1/6, 43/180, 3/10, 197/560.
5. Identity (I): NOT FOUND stated anywhere. Nearest hit: Finch (arXiv:1601.04937 p.14) notes
   the Gaussian instance as a curiosity ("It is interesting that sqrt(3) is twice the expected
   area of a planar Gaussian triangle") — and that is about LEBESGUE area, a logically
   different statement from (I), which is about the mu-content. Both verified numerically here.
   Verdict: treat (I) as folklore-or-known; DO NOT CLAIM AS NEW.
6. (I) is distribution-free: verified by MC for a standard Gaussian, a skewed anisotropic
   Gaussian mixture and a heavy-tailed t_3 law (|z| < 1.7) — no convexity, no uniformity.
   `results/A4_gaussian_identity_check.json`.

## 2026-08-19 00:15 — [server] BUG FOUND AND FIXED in the landscape scan (not a counterexample)
The first n=5 landscape run reported `max over <= 3-gons: F = 0.076256 *** EXCEEDS TRIANGLE
BOUND ***`. Since every triangle is an affine image of every other and F is affine invariant,
F must be the CONSTANT 5/72 on 3-gons — so this was a bug by construction, and it was:
Nelder-Mead had wandered onto an almost-collinear triangle (area 2.4e-5, edge lengths 4.1,
1.4, 2.7) where `cov_det_area`'s float64 `E[xx^T] - mu mu^T` cancellation destroyed det Sigma
(E[A_3^2] came out 0.00708 instead of 1/72 = 0.01389). Route P itself was fine there
(E[A_3] = 0.08333333332934, T_0 - 2 = 8.6e-12).
Fixes in `polygon_exact.py` / `n5_landscape.py`:
 * `cov_det_area_exact` — the covariance determinant in exact `Fraction` arithmetic from the
   float vertices (no cancellation at all); on the offending triangle it returns 1/72 exactly.
 * `whiten` — put every polygon in isotropic position (centroid 0, Sigma = I) before
   evaluating anything. Legitimate because every quantity here is affine invariant, and it
   makes route P well-conditioned too.
 * route P's own checks `|T_0 - 2| < 1e-9`, `|E[N_3] - 3| < 1e-9` are now enforced per body;
   failure returns NaN instead of a number.
 * a self-test that runs BEFORE the scan: F over 200 random affine images of a triangle and
   200 of a parallelogram (condition numbers up to ~1e6) must match 5/72 and (1-49/144)/10.
After whitening, the offending triangle gives F = 0.069444444444444 = 5/72 exactly.
LESSON (for LESSONS_LEARNED): an affine-invariant quantity is its own regression test —
evaluate it on random affine images and require the invariance to machine precision.

## 2026-08-19 00:19 — [server] B DONE: exact P_5 for regular m-gons + landscape verification
`docs/N5_LANDSCAPE.md`, `results/n5_landscape.json`, `results/B1_mgon_exact_table.json`.
CLOSED FORM (Alikoski 1939 + identity (II)), verified three ways:
    P_5(regular m-gon) = 1 - 5 (15 cos^2 w + 92 cos w + 76) / (36 m^2 sin^2 w),  w = 2 pi/m
  m=3 -> 11/36 (Valtr), m=4 -> 49/144 (Valtr), m->oo -> 1 - 305/(48 pi^2) (Marckert).
  m=5: 7/12 - 47 sqrt5/450 = 0.349788455683355
  m=6: 1373/3888                = 0.353137860082305
  m=8: 1469/2304 - 115 sqrt2/576 = 0.355235139456761
  m=10: 461/720 - 229 sqrt5/1800 = 0.355800240640305
  m=12: 3439/5184 - 115 sqrt3/648 = 0.356001785693795
LANDSCAPE (561 bodies, all exact via route P, affine-whitened): regular m-gons m=3..40,
truncated triangles, circular caps, stadiums, triangle->disk Minkowski interpolations,
pushed-vertex polygons, 400 random convex polygons, half-disk.  **0 violations** of
61/(96 pi^2) <= F <= 5/72 at tolerance 1e-12. Direct Nelder-Mead optimisation over convex
polygons with <= m vertices returns the triangle value 5/72 as the maximum for every
m in {3,4,5,6,8,10}, and the regular m-gon value as the minimum for m = 4,5,6,8 exactly
(m=10 under-converged: 0.0644406 vs the decagon's 0.0644200 — a limitation, logged as such).
Self-test before the scan: F on 200 random affine triangles and 200 random parallelograms
reproduces 5/72 and (1-49/144)/10 to 5.3e-16 / 7.2e-16.
CAVEAT, prominent in the doc: this is a VERIFICATION of Marckert & Rahmani's 2021 theorem,
not a test of an open conjecture. Only n >= 6 is still open.

## 2026-08-19 00:19 — [server] C: P_6(disk) fully assembled and cross-checked
`src/disk_EA5_exact.py`, `src/C_disk_n6_moments.py`.
NEW (this session, sympy, from E[N_6] = 30 E[c^4] and the master formula):
    E[N_6](disk)  = 6 - 175/(12 pi^2) + 23023/(1152 pi^4)
    E[A_5](disk)  = 7 (2400 pi^2 - 3289) / (6912 pi^4) = 0.2120720740367009709...
  Independent confirmation: route P on polygonal approximations + Richardson gives
  0.212072074036701 -- agreement 2.9e-17.
FROM MARCKERT 2017's P^D_{6,m} (read verbatim) plus the classical Efron-Buchta identities:
    E[A_3^3](disk) = P_{6,3}/20             = 1001/(6400 pi^4)
    E[A_4^2](disk) = P_{6,4}/15 + 4E[A_3^3] = (2400 pi^2 + 31031)/(19200 pi^4)
CONSISTENCY (all exact, all zero):
  * Marckert's P_{5,m} and P_{6,m} each sum to 1 exactly (so the PDF transcription is right).
  * our identities give P_{5,3} = 10E[A_3^2] = 15/(16 pi^2) and P_{5,4} = 5E[A_4]-20E[A_3^2]
    = 65/(12 pi^2) -- exactly Marckert's values.
  * E[N_6] computed from Marckert's P_{6,m} equals our sympy c^4 value EXACTLY (diff 0),
    two completely different derivations meeting.
  * 1 - 6E[A_5] + 15E[A_4^2] - 20E[A_3^3] = 1 - (146400 pi^2 - 473473)/(11520 pi^4) exactly,
    i.e. Marckert's P_6.  A 2e8-sample MC check of the three moments is running.

## 2026-08-19 00:24 — [server] NEGATIVE RESULT: no Alikoski-shaped closed form for E[A_5](m-gon)
`src/mgon_EA5_fit.py`, `results/C_mgon_EA5_fit.json`. Alikoski's four-point answer has the
shape (quadratic in C)/(36 m^2 S^2) with C = cos w, S = sin w, w = 2pi/m, and the disk limit
of E[A_5] has both a 1/pi^2 and a 1/pi^4 part, so the natural ansatz is
    E[A_5] = a_0 + P(C)/(m^2 S^2) + Q(C)/(m^4 S^4),   deg P <= 3, deg Q <= 4   (10 unknowns).
Solved exactly at 60 dps on m = 3..12, then tested on m = 13,14,15,16,20,24,30,40:
residuals 2.7e-12 ... 6.6e-10 and GROWING with m (they would be ~1e-45 if the ansatz held);
fitted coefficients blow up to ~1e6, i.e. the basis is nearly degenerate on the fit points.
VERDICT: E[A_5](regular m-gon) is NOT of this shape. Reported as a failure, not massaged.
(High-precision E[A_5] values for m = 3..60 are in `results/C_mgon_EA5.json`; they do
reproduce Buchta's 43/180 (m=3) and 79/360 (m=4) exactly, and satisfy identity (I')
E[A_6] = 3E[A_5] - 5E[A_3] to 1e-48 for every m.)

## 2026-08-19 00:31 — [server] identity (II) verified for a NON-UNIFORM law (Gaussian)
`src/gaussian_P5_test.py`, `results/A3_gaussian_P5.json`, 1.2e8 samples. For a general
absolutely continuous mu the two moments are plain containment frequencies:
E[A_3] = P(X_4 in conv(X_1,X_2,X_3)) and E[A_3^2] = P(X_4 AND X_5 in conv(X_1,X_2,X_3)),
so no area computation enters at all. Standard bivariate Gaussian:
  E[A_3]   = 0.087741100 +- 2.4e-05
  E[A_3^2] = 0.014652233 +- 9.7e-06
  P_4 direct MC = 0.649020550 +- 3.3e-05
      vs 1 - 4E[A_3] (Sylvester/Efron)              z = -0.16
      vs the classical 6/pi * arcsin(1/3) = 0.649040688   z = -0.60   [classical anchor OK]
  P_5 direct MC = 0.269051383 +- 4.5e-05
      vs identity (II) 1 - 10E[A_3] + 10E[A_3^2] = 0.269111333   z = -0.26
So (II) is not an artefact of uniformity or of bounded convex support. (The covariance form
P_5 = (5/2)P_4 - 3/2 + 15 detSigma/|K|^2 is uniform-specific -- Lemma 2 needs A_3 = |det|/(2|K|)
-- but the moment form 1 - 10E[A_3(1-A_3)] is general.)
NOTE: P_5(standard Gaussian) = 0.26905 +- 0.00005 is, as far as this session's search went,
not a value we found tabulated anywhere. It is exactly 1 - 10E[A_3] + 10E[A_3^2] with
E[A_3] = (1 - (6/pi) arcsin(1/3))/4 exact; only E[A_3^2] is missing a closed form. UNVERIFIED
as a novelty claim -- not searched for specifically.

## 2026-08-19 00:32 — [server] C VERIFIED: all six disk moments confirmed by Monte Carlo
`results/C_disk_n6_moments.json`, 2e8 samples per k, Jarvis-march hull areas (route M code,
structurally independent of both the sympy derivation and Marckert's paper):
  E[A_3]   MC 0.073879359 +- 4.1e-06  exact 35/(48 pi^2)                    z = -0.16
  E[A_3^2] MC 0.009498935 +- 9.5e-07  exact 3/(32 pi^2)                     z = +0.08
  E[A_3^3] MC 0.001605706 +- 2.4e-07  exact 1001/(6400 pi^4)   [Marckert]   z = +0.18
  E[A_4]   MC 0.147760203 +- 7.6e-06  exact 35/(24 pi^2)                    z = +0.02
  E[A_4^2] MC 0.029255291 +- 2.9e-06  exact (2400 pi^2+31031)/(19200 pi^4)  z = -0.60
  E[A_5]   MC 0.212074089 +- 5.6e-06  exact 7(2400 pi^2-3289)/(6912 pi^4)   z = +0.36
  P_6 reassembled from the MC moments: 0.134270702 +- 5.2e-05 vs Marckert 0.134309386, z = -0.74
Every extracted moment stands up. Section C is closed: P_6(disk) is known (Marckert), and all
three of its ingredients now have exact values that three independent routes agree on.

## 2026-08-19 00:33 — [server] IMPORTANT: docs/PROBLEM_B.md was briefing a SOLVED problem
Found `docs/PROBLEM_B.md` (dated 2026-08-18 23:45, not written by this session — it appeared
in `docs/` after the session started, presumably synced by the operator for another agent).
It poses the n=5 extremal conjecture as an open problem ("To our knowledge n = 5 is open")
and asks for a proof or disproof of exactly the two inequalities that Marckert & Rahmani
proved in 2021. Per TASK.md ("If a claim in this brief turns out wrong, say so in PROGRESS.md
and fix it everywhere"), I have PREPENDED a clearly-marked correction banner to that file
rather than editing its body: it flags the theorem, the Marckert-Morin verbatim confirmation,
the fact that P_5(disk) is in Marckert 2017, the exact pentagon/hexagon values replacing its
Monte-Carlo estimates, and one misleading sentence about the direction of the two extremal
results. The original text is left intact below the banner.
ANYONE READING PROBLEM_B.md SHOULD READ THE BANNER FIRST.

## 2026-08-19 00:38 — [server] NEW: a 3-D analogue of identity (I), proved and tested
Arose from checking a speculative sentence before writing it into VERDICT.md — the probe
(`src/dim3_probe.py`) showed E[A_5]/E[A_4] ~ 2.5 for the 3-D ball, cube and simplex, and
there is a reason. The planar proof of (I) uses #vertices = #edges, false in R^3; but for a
SIMPLICIAL 3-polytope Euler V - E + F = 2 with 3F = 2E gives V = 2 + F/2, and generic points
give a simplicial hull. With the R^3 Renyi-Sulanke facet count
E[F_n] = C(n,3) E[c^{n-3} + (1-c)^{n-3}] = 2 C(n,3) E[c^{n-3}] (c = mass on one side of the
plane through three points; c =_d 1-c by the same transposition argument), Efron gives
    E[A_{n-1}] = 1 - 2/n - C(n,3) m_{n-3}/n,
    E[A_3] = 0 (correct), E[A_4] = 3/5 - 2 m_2, E[A_5] = 2/3 - (10/3) m_3 = 3/2 - 5 m_2.
    ==>  E[A_5] = (5/2) E[A_4]   in R^3, for every absolutely continuous law.   (I_3)
    and  E[A_7] = (7/2) E[A_6] - (35/4) E[A_4].
CONSEQUENCES (five-point hull volumes free from the four-point ones):
    E[A_5](ball) = (5/2)(9/715) = 9/286 = 0.0314685314...
    E[A_5](cube) = (5/2)(3977/216000 - pi^2/2160) = 3977/86400 - pi^2/864 = 0.0346069394...
TEST (`src/dim3_identity.py`, 3e7 samples/body; N_5 and N_6 from exact orientation
predicates only -- a point is a non-vertex iff it is in conv of four of the others
(Caratheodory), so no hull construction and no floating volume enters):
  ball    E[A_4]=0.012600020+-1.1e-05 (anchor 9/715, z=+1.19); E[A_5]=0.031468850+-1.2e-05
          (vs 9/286: z=+0.03);  (I_3) residual -3.12e-05 +- 2.4e-05, z=-1.30
  cube    E[A_4]=0.013844173+-1.0e-05 (anchor Zinani, z=+0.14); E[A_5]=0.034581361+-1.5e-05
          (vs 3977/86400 - pi^2/864: z=-1.71);  (I_3) residual -2.91e-05 +- 2.0e-05, z=-1.49
  simplex E[A_4]=0.017390713+-1.2e-05; E[A_5]=0.043506944+-1.3e-05;
          (I_3) residual +3.02e-05 +- 2.5e-05, z=+1.21
All |z| < 1.5.  NOTE: the first, cruder probe using scipy ConvexHull volumes suggested a
2.5-sigma excess for the ball; that was an artefact of the hull-volume estimator plus a
1.5-sigma low E[A_4], and it disappeared with the exact-predicate estimator. The lesson is
the same one as at 00:15 — prefer estimators that are exact combinatorics over ones that
build floating-point geometry.
NOVELTY: NOT FOUND stated (searches in LITERATURE.md); same caveat as (I) -- it is one step
from classical ingredients, so do not claim as new without a proper search.

## 2026-08-19 01:09 — [server] A3 COMPLETE: identities (I) and (II) on nine bodies
`results/A3_identity_mc.json`. Reference E[A_3] exact (square, triangle), exact route P
(pentagon, hexagon, octagon) or route-P Richardson (disk, 3:1 ellipse, half-disk, stadium);
E[A_3^2] exact in every case. Direct tester 2e8 samples for P_4 and P_5, moments 1e8.
  body       ref                       z(I)   z(P_4)  z(II)   testers agree
  square     exact rational           -0.06   +1.51   +0.40   yes
  triangle   exact rational           +1.68   -0.32   +2.09   yes
  disk       route P Richardson       +1.25   -0.44   +0.47   yes
  pentagon   route P exact polygon    +1.32   +0.65   -1.48   yes
  hexagon    route P exact polygon    +1.06   +1.71   -0.14   yes
  octagon    route P exact polygon    +1.16   +0.03   +0.96   yes
  ellipse3   route P Richardson       +1.25   -0.44   +0.47   yes
  halfdisk   route P Richardson       -0.11   +0.93   +0.12   yes
  stadium    route P Richardson       -1.04   +1.33   +2.47   yes
max |z| = 2.47 over 27 tests; 0 above 3; testers A and B agree sample-for-sample everywhere.
NEW high-precision values for bodies with no literature value we found:
  half-disk : E[A_3] = 0.076512497523552, P_4 = 0.693950009905794, P_5 = 0.341069936631669
  stadium(a=1): E[A_3] = 0.074403975610728, P_4 = 0.702384097557086, P_5 = 0.352809265308243
  (note the half-disk has P_4 BELOW the square's 25/36 = 0.694444 but P_5 ABOVE 49/144.)
TWO CAVEATS, recorded rather than glossed:
 * `ellipse3` is NOT an independent test. sample_ellipse3 is the disk sampler scaled by
   diag(3,1) with the same seed, so it draws the affine image of the SAME point sets and
   reproduces the disk's z-scores digit for digit. That is a valid (and exactly passing)
   affine-invariance test of the whole pipeline, but it is one body, not two.
 * The 27 z-scores have mean +0.62 (sd 0.96), which is 3.2 sigma from 0 if they were
   independent. They are not: every body used the SAME seeds (202 for n=4, 303 for n=5), so
   different bodies consume one uniform stream through different transformations and their
   fluctuations are correlated. A dedicated check is running (`src/bias_check.py`):
   (A) square P_4 against the exact 25/36 with 12 INDEPENDENT seeds -- the real bias test;
   (B) the whole body set re-run with one fresh common seed, to see whether the common
   offset simply moves. Result in the next entry.

## 2026-08-19 01:10 — [server] the mean-z tilt was a SEED ARTEFACT, not a bias (settled)
`src/bias_check.py`, `results/A3_bias_check.json`.
(A) THE ACTUAL BIAS TEST — square P_4 against the exact 25/36, 1e8 samples x 12 INDEPENDENT
    seeds: z = -0.14, +0.47, +1.89, -0.02, -0.21, +0.09, +1.38, -0.91, +0.29, -1.84, -1.19
    (+ seed 9001), mean z = -0.106 +- 0.289, sd = 1.065.
    => no evidence of bias, and the batch/binomial standard errors are correctly sized
       (an sd of 1.065 on 12 draws is exactly what N(0,1) should give).
(B) THE CORRELATION TEST — all eight bodies, n=5, 1e8 samples, one common seed:
       old seed 303:   +0.91 +0.22 -0.52 -1.30 +0.14 +1.18 -0.15 +2.08   mean +0.317
       fresh seed 555: -0.42 -0.90 +1.86 +0.72 -0.08 -0.56 -0.33 -1.84   mean -0.193
    The common offset MOVES with the seed and changes sign, exactly as predicted: with a
    shared seed the bodies consume one uniform stream through different transformations, so
    their z-scores are positively correlated and the ensemble mean is not a bias statistic.
    The stadium's +2.47 in the main run is a fluctuation of that kind -- it goes to -1.84 on
    a fresh seed.
CONCLUSION: the 27 z-scores of the 01:09 entry stand as reported; their mean of +0.62 carries
no information about bias and should not be read as evidence of one. LESSON worth exporting
to LESSONS_LEARNED: when sweeping many configurations, VARY THE SEED PER CONFIGURATION,
otherwise the z-scores are correlated and the ensemble looks either falsely biased or falsely
clean. Test bias with many seeds on ONE configuration with a known exact answer.

## 2026-08-19 01:17 — [server] VERDICT.md written; session complete, stopping as instructed
All of TASK.md sections A, B and C are done and written up. No server jobs from this campaign
remain running. Deliverables: THEOREMS.md, LITERATURE.md, docs/N5_LANDSCAPE.md, VERDICT.md,
a correction banner on docs/PROBLEM_B.md, 38 JSONs in results/, and 11 new src modules.
Per the brief ("When A, B are done ... write VERDICT.md, then STOP. Do not start new research
directions"), stopping here. The single highest-value follow-up is to obtain Marckert &
Rahmani, Mathematika 67 (2021) 860-884 — it decides the novelty status of identity (II) and
of the regular-m-gon P_5 closed form, and it is the only source this sweep could not reach.

## 2026-08-19 01:18 — [server] one redundant job killed, its data salvaged
`mgon_EA5.py` was still running after 1h, stuck in its trailing `mpmath.identify()` calls on
algebraic numbers (a slow, open-ended constant-recognition attempt). Its useful output — the
50-digit E[A_5] table for m = 3..60 and the identity-(I') residuals — was already complete in
`results/C_mgon_EA5.log` and has been parsed into `results/C_mgon_EA5.json` (19 rows). The
closed-form question it was trying to answer was settled decisively and separately by
`mgon_EA5_fit.py` (00:24 entry: the Alikoski-shaped ansatz FAILS). Killed it rather than let
it run unattended. No campaign processes remain.


## 2026-08-19 06:14 — local session back; overnight VERDICT read; external n=5 proof integrated
Overnight session (23:40-01:17) finished: VERDICT.md, THEOREMS.md, LITERATURE.md, N5_LANDSCAPE.md,
route P (exact polygon Sylvester quantities), P_5(regular m-gon) closed form, disk n=6 assembled
from Marckert 2017, 3-D analogue E[A_5] = (5/2)E[A_4]. Headline: P_5(disk) is in Marckert 2017
Table (7); the n=5 extremal conjecture is Marckert-Rahmani 2021 (only n>=6 open).
06:20  Operator relayed an external model's proof of Problem B (section lemma). Verified
       (src/section_lemma_check.py): plateau + quadratic law on a fibre; Steiner deficit on the
       half-disk and shaking gap on the disk both = 0.001511838 +- 5.5e-11 vs target
       F(half-disk)-F(disk) = 0.001511837564 from route P + exact covariance. PASS.
06:30  docs/N5_PROOF.md written; pointers added to THEOREMS.md (sec 8), N5_LANDSCAPE.md (sec 6),
       PROBLEM_B.md, LITERATURE.md (sec 7), VERDICT.md (addendum). Server mirror pulled into
       the local repo (78 result files). Local literature agents from 23:25 produced no docs/lit_*.md
       (their reports never landed) -- superseded by the server session's LITERATURE.md.

06:35  Marckert's paper located on arXiv (1511.03658) via the Semantic Scholar DOI record; full text
       extracted (docs/marckert_rahmani_1511.03658.txt) and read. Identity (II) / E[A_4]=2E[A_3] /
       covariance form / regular polygons: NOT present. Its n=5 proof: comb recursion + polynomial
       comparison in the symmetry defect (quadratic for n=5; degree 2 and 4 for n=6, unresolved).
       LITERATURE.md sec 8, N5_PROOF.md sec 6, VERDICT.md addendum updated. Journal version unread.


## 2026-08-19 (morning cont.) — n = 6 ingredient table for triangle/square/disk (EXACT)
Decomposition P_6 = 1 - 6E[A_5] + 15E[A_4^2] - 20E[A_3^3] validated vs Marckert's disk value
(n6_bp_moments.py). New exact tool: a Blaschke-Petkantschin WIDTH-FUNCTION integral giving
E[A_3^k] for any convex body (n6_bp_polygon.py; disk k=1,2,3 to >30 digits; triangle/square
k=1,2 exact). Results (PSLQ from 40-digit integrals, MC-confirmed rel ~1e-4):
  E[A_3^3]: triangle 31/9000, square 137/72000, disk 1001/(6400 pi^4).
Same machine (via E[c^j]) gives E[A_5] = 43/180 (triangle), 79/360 (square) [Buchta 1984] and
re-confirms E[A_4]=2E[A_3]. Backing E[A_4^2] out of the decomposition + Valtr and CONFIRMING by
MC (mc_a4sq.py, corrected hull-4 area = max over 3 quad-orders + 4 drop-one triangles):
  E[A_4^2]: triangle 181/4500, square 859/27000, disk (2400 pi^2+31031)/(19200 pi^4).
All three P_6 loops close (n6_results.py: ALL LOOPS OK). Full write-up docs/N6_INGREDIENTS.md,
table results/n6_ingredients.json. Bugs caught: (a) phi-subdivision must include DIAGONAL
directions not just edge normals (square failed at high precision until fixed); (b) hull-area
of 4 points needs the drop-one triangles for interior-point configs (first MC was 5-7% low).
OPEN: an INDEPENDENT exact route to E[A_4^2 & convex] for a general body (via hull area =
(1/2)|diag1 x diag2|) -> would give first-principles P_6 for regular polygons.


## 2026-08-19 08:10 — server session opens TASK_N6.md (two-chord polygon integral)
08:09  Read README / docs/N6_INGREDIENTS.md / TASK_N6.md. Re-ran the disk reference
       src/n6_twochord_disk.py: E[A_4^2 & convex](disk) = 0.02283435, rel.diff 1.53e-12 vs
       (2400 pi^2+19019)/(19200 pi^4). Machinery + normalisation (1/V^6) confirmed.
08:10  Plan for the polygon version: change variables from (p1,p2) to the crossing point X
       (dp1 dp2 = |sin(phi1-phi2)| dX), giving
         I = int_0^pi int_0^pi sin^2(D)|sin(D)| F(phi1,phi2) dphi1 dphi2,
         F = int_K G1(X) G2(X) dX,  G_i(X) = ((u_i+v_i)^5-u_i^5-v_i^5)/10,
       u_i,v_i = the two ray distances from X to dK along +-t(phi_i).  Inside a cell
       (p1-strip of direction 1) x (p2-strip of direction 2) each of u1,v1,u2,v2 is LINEAR in X,
       so G1*G2 is a degree-10 polynomial and the inner 2-D integral is done EXACTLY by a
       degree-11 Duffy-Gauss rule on a fan triangulation of the clipped cell.  Only the 2-D
       angular integral is numerical, panelled at every critical angle (vertex-pair directions
       and their perpendiculars) plus phi2 = phi1 (the |sin|^3 kink).  Integrand is everywhere
       >= 0, so no cancellation.
08:25  src/n6_twochord_polygon.py written and VALIDATED. Independently of Valtr:
         E[A_4^2 & convex](triangle) = 0.026444444444  vs 119/4500   rel 1.7e-15  (nphi=24)
         E[A_4^2 & convex](square)   = 0.024203703704  vs 1307/54000 rel 6.6e-13  (32x24)
       Cross-checks: (a) the raw right-triangle / unit-square vertex lists with NO symmetry
       reduction (q=1) give the same values as the regular 3-gon (q=3) / 4-gon (q=2) --
       validates the rotation-symmetry fundamental-domain factor; (b) the answer is exactly
       invariant under raising the Duffy-Gauss order ng from 6 to 8 to 10, confirming the inner
       2-D integral really is exact (degree-10 polynomial per cell), so the ONLY error is the
       2-D angular quadrature. STEP 3 OF TASK_N6 PASSED.
08:55  First-principles P_6(regular m-gon) table computed, m = 3..12 (src/n6_mgon_P6.py,
       results/n6_twochord_polygon.json).  ANCHORS: m=3 -> 0.10111111111111069 (91/900,
       err 4e-16), m=4 -> 0.1224999999999184 (49/400, err 8e-14); P_6 increases monotonically
       in m towards the disk 0.134309386.  Ingredients PSLQ'd exactly in Q(cos 2pi/m):
         E[A_3^3]: 31/9000 (m=3), 137/72000 (4), 32/28125+19sqrt5/75000 (5), 57709/34992000 (6),
                   193/230400+53sqrt2/96000 (8), ... , 6103/7464960+1183sqrt3/2592000 (12)
         E[A_5]  : 43/180, 79/360, 221/1500+34sqrt5/1125, 149347/699840,
                   4531/36864+5837sqrt2/92160, ..., 253843/2239488+79843sqrt3/1399680
       (E[A_3] = (9c^2+52c+44)/(36 m^2 s^2) [Alikoski] and E[A_3^2] = (2+c)^2/(24 m^2 s^2)
       re-verified to 30 digits by the same width route -- framework check.)
09:20  TWO NUMERICS FIXES that mattered:
       (a) G = ((u+v)^5-u^5-v^5)/10 = u v (u+v)(u^2+uv+v^2)/2 -- the factored form is
           cancellation-free; in sliver cells u/v reaches 1e10 and the naive difference lost
           ~10 digits.
       (b) graded angular panels (geometric refinement towards every critical angle, where
           t(phi) becomes parallel to an edge and F is only finitely smooth).
       Plus an 80-bit (np.longdouble) arithmetic path with mpmath-generated Gauss nodes.
       Square now reproduces 1307/54000 to 8e-16 rel (was 8e-8 at the first grid).
09:35  E[A_4^2 & convex] to ~2e-17 absolute for m=3..12 -> exact values by a SMOOTH-DENOMINATOR
       scan (test D = 2^a 3^b 5^c 7^d, demand D*x integral / in Z[c]).  Calibration: recovers
       119/4500 (m=3) and 1307/54000 (m=4).  NEW: E[A_4^2 & convex](hexagon) = 403891/17496000,
       hence  E[A_4^2](hexagon) = 57701/1944000  and  P_6(regular hexagon) = 461299/3499200
       = 0.1318298468221308.  (3499200 = 2^6 3^7 5^2, same family as the other hexagon
       denominators 34992000 = 2^7 3^7 5^3 and 699840 = 2^6 3^7 5.)
09:30  CORRECTION + three precision bugs found and fixed (the 08:35 hexagon value was wrong):
       (i)   `cossin` rebuilt the angle as mp.mpf(repr(float(phi))) -- a decimal round-trip that
             re-rounds, costing 1.5e-17;  mp.mpf(python_float) is exact, use that.
       (ii)  the phi integration limits (pi, pi/q) and the critical angles were float64, so the
             OUTER limits of the integral carried a 1e-16 relative error.
       (iii) the decisive one: inside the cell loop, `affine()` did
             `den = float(nu[jj] @ d); return (float(h[jj])/den, ...)` -- two float() casts that
             computed the affine coefficients of u,v in DOUBLE precision inside the 80-bit path,
             biasing every F by ~8e-16 relative.  Diagnosed by writing an exact mpmath reference
             for a single F(phi1,phi2) (same cell algorithm, mp.mpf throughout) and comparing:
             ld/mpmath relative difference was a consistent +7..9e-16 before, ~1e-19 after.
       After the fixes: E[A_4^2 & convex] reproduces 119/4500 to 3.2e-19 relative and
       1307/54000 to 1.3e-17 relative -- i.e. ~19 significant digits, the longdouble limit.
       LESSON (again): a plateau that does not move with grid refinement OR with the arithmetic
       precision is a bug, not a quadrature limit; and the way to find it is an
       independent-precision reference for the innermost quantity, not more refinement.
       Re-run with the fixes: E[A_4^2 & convex](hexagon) = 0.023084762231367169649, and
       403891/17496000 = 0.023084762231367169639 -- so the 08:35 hexagon candidate SURVIVES,
       now confirmed to 1e-20 instead of resting inside a 2e-17 error bar.  The m=5 candidate
       from the same scan also survives (now matched to 1e-21).

## 2026-08-19 10:15 — TASK_N6 done: P_6(regular m-gon) from first principles
10:15  Final table (src/n6_mgon_final.py -> results/n6_mgon_P6_final.json;
       E[A_4^2 & convex] from results/n6_tc_hiprec.json, three independent (nphi,grade)
       settings in 80-bit agreeing to ~1e-20; E[A_5], E[A_3^3] from the width route):
         m :        P_6                       exact
         3 : 0.101111111111111111239   91/900                              [Valtr -- OUTPUT now]
         4 : 0.122499999999999986519   49/400                              [Valtr -- OUTPUT now]
         5 : 0.129248382075802995501   8941/22500 - 1349 sqrt5/11250            *** NEW ***
         6 : 0.131829846822130772901   461299/3499200                           *** NEW ***
         7 : 0.132961474013140063212   (deg-3 field, not identified)
         8 : 0.133516325200578055307   30103/61440 - 116141 sqrt2/460800        *** NEW ***
         9 : 0.133813207438104335867   (deg-3 field)
        10 : 0.133983396871405807140   (deg-2 field, denominators too large for 19 digits)
        11 : 0.134086526178116899361   (deg-5 field)
        12 : 0.134151930266256898067   (deg-2 field, ditto)
       inf : 0.134309386357109938800   1 - (146400 pi^2 - 473473)/(11520 pi^4)  [Marckert]
       Exact-vs-numeric residuals for m=3,4,5,6,8: 1.3e-19, 1.3e-17, 1.4e-20, 1.5e-19, 1.1e-19.
       New E[A_4^2 & convex]: 1769/112500 + 577 sqrt5/168750 (pentagon), 403891/17496000
       (hexagon), 35743/2764800 + 48793 sqrt2/6912000 (octagon).
10:20  `python n6_twochord_polygon.py --validate` (16 checks: eqtri/right-triangle/diamond/unit
       square x {f64, ld, ld nphi=48, ld ng=10}) -> ALL ANCHORS OK, worst rel 9.2e-16 (f64),
       best 2.1e-20.  results/n6_twochord_polygon.json.
10:25  NEGATIVE RESULT 1: no P_5-style closed form.  Fitting X = A(c)/S + B(c)/S^2, A,B in Q[c],
       S = m^2 sin^2 w, to 11 points (m=3..12 + disk) each good to ~1e-18: best 9-parameter
       maxrelres 8e-8, decaying smoothly with parameter count -> interpolation, not discovery.
       A wider automated search over sum_j a_j c^j/(m^a sin^b w) (a,b<=8, <=3 groups, <=10
       parameters) found nothing for E[A_3^3] either.  Consistent with the overnight session's
       independent failure on E[A_5](m-gon) (results/C_mgon_EA5_fit.json, "fails").
       NEGATIVE RESULT 2: the leading large-m correction is 1/m^4, not 1/m^2.  For n=5 the exact
       closed form gives P_5(disk)-P_5(m) = (7 pi^2/18)/m^4 + O(1/m^6) (the 1/m^2 term cancels
       identically; 7pi^2/18 = 3.8381795, confirmed numerically).  See 11:40 for n=6.
10:30  Extremal-conjecture data point: P_6 strictly increasing in m from 91/900 (triangle) to the
       disk 0.1343093863571 -- the ordering the OPEN n>=6 conjecture predicts, now resting on
       EXACT values at m=3,4,5,6,8 rather than Monte Carlo.
10:35  Independent cross-check: direct convex-position Monte Carlo (src/mc_mgon_n6.py, the numba
       hull tester -- no shared code with the two-chord / width pipeline), 6e8 samples:
       pentagon 0.12924843 +- 1.37e-05 vs 0.12924838 assembled, z = +0.00;
       hexagon  0.13184907 +- 1.38e-05 vs 0.13182985, z = +1.39;
       octagon  0.13353532 +- 1.39e-05 vs 0.13351633, z = +1.37.  (The two +1.4's share seed 7,
       so they are one correlated fluctuation, not two.)  results/n6_mgon_mc_check.json.
10:40  Write-ups: docs/N6_LANDSCAPE.md (new, sibling of N5_LANDSCAPE.md), docs/N6_INGREDIENTS.md
       ("RESOLVED" section), VERDICT.md addendum, README.md table + contents.
11:40  Large-m extension (src/extend_a4.py; m = 14, 16, 18, 20; 19-45 min each, the mpmath
       width route being the bottleneck) -> P_6 = 0.1342243249258609, 0.1342595003080803,
       0.1342782328411357, 0.1342889420922087.  Richardson fits in 1/m^2 with 5-10 terms now
       agree to 1.3e-9:  a_4 = 3.2738126823 +- 1.3e-9  (a_4/pi^2 = 0.3317065760).
       AND A CANDIDATE DIED: at 6 digits (m <= 12) PSLQ over {1, pi^2, pi^-2} had returned the
       very clean  a_4 = (4 pi^2 + 218 + 76/pi^2)/81  -- one denominator, three small numerators,
       matching to 8e-10.  The sharpened a_4 is 3.6e-9 away, i.e. 2.9 spreads: REFUTED.  The
       sharpened value in turn produces a different "clean" relation at every basis and height
       (11/36, 19/180, 271/180 at height 300; 4/7, -57/14, 247/14, -8 with a pi^-4 term).  Nine
       digits is not enough for a 3-4 term relation.  a_4 stays UNIDENTIFIED -- exactly the
       failure mode the campaign's "PSLQ only from many digits" rule exists to catch.
       results/n6_mgon_large_m.json; table rows merged into results/n6_mgon_P6_final.json.
11:50  Final state: no campaign processes left running.  Deliverables of TASK_N6:
         src/n6_twochord_polygon.py   (--validate: 16 anchor checks, ALL OK)
         src/n6_mgon_final.py, src/pslq_n6.py, src/mc_mgon_n6.py, src/extend_a4.py,
         src/tc_hiprec.py, src/tc_sweep.py, src/search_form.py, src/fit_n6_mgon.py
         results/n6_twochord_polygon.json  (anchor validation + per-m values with error bars)
         results/n6_mgon_P6_final.json     (the assembled P_6 table, m = 3..20 + disk)
         results/n6_tc_hiprec.json, n6_tc_sweep.json, n6_mgon_width_ingredients.json,
         results/n6_mgon_mc_check.json, n6_mgon_large_m*.json
         docs/N6_LANDSCAPE.md (new), docs/N6_INGREDIENTS.md ("RESOLVED"), VERDICT.md addendum,
         README.md (anchor table + contents).

## 2026-08-19 11:46 (local) — pulled server n6 deliverables, independently re-verified
Server TASK_N6 completed ~11:50 (ran ~3h35m from 08:05). Pulled into local repo. Independent
cross-check with the LOCAL width route (src/n6_bp_polygon.py): E[A_3^3](hexagon) and E[A_5](hexagon)
match the server's to ~1e-27; assembling P_6 from MY E[A5],E[A3^3] + the SERVER's two-chord
E[A_4^2&convex] gives exactly 461299/3499200 (match 1.5e-19). Two-chord polygon route confirmed.
Server session left an UNSUBMITTED line in tmux: "push the m=10 and m=12 identification with more
digits" (m=10 in Q[sqrt5], m=12 in Q[sqrt3] have 19-20 digit decimals but not yet exact closed
forms) -- a reasonable follow-up, not yet run.

## 2026-08-19 12:31 (local) — handed n>=6 extremal scan + R^d to server (TASK_N6_EXTREMAL.md)
Local exact P_6-for-any-polygon driver src/n6_p6_body.py validated (right triangle 91/900 to 1e-15)
but ~15-60s/body -> too slow for a broad shape scan. The conjecture window [91/900, 0.134309] is
wide, so direct convex-position MC (4-5 digits) is the right scan tool; handed to the server.

## 2026-08-19 (server) — TASK_N6_EXTREMAL: the n >= 6 extremal scan

12:44  Part 1 step 1 DONE.  `src/mcp6.py`: fast direct convex-position MC of P_6 AND P_7 for an
       arbitrary convex polygon, numba-parallel, 2e8 samples in 4.1 s on 16 threads (~5e7
       samples/s; each sample draws 7 points, the first 6 giving P_6 and all 7 giving P_7 in one
       pass).  Sampling = centroid fan triangulation + area-weighted binary search; convex
       position tested by "sort by angle about the sample centroid, then require a left turn at
       every vertex" (equivalent to hull-size == n; angles compared via the transcendental-free
       monotone diamond pseudo-angle, turns by exact float orientation determinants).  RNG is
       counter-based splitmix64 re-seeded per chunk, so output depends only on
       (seed, nchunk, nsamp) and NOT on the thread count.
12:44  VALIDATION (`mcp6.py --validate --samples 4e8`, results/n6_mc_validation.json): 24 checks,
       ALL OK.  17 exact anchors, every z in [-1.4, +1.1]:
         P_6: 91/900 (right / equilateral / 1:100-sliver triangle), 49/400 (unit square and a
              sheared parallelogram), the exact regular m-gon table m = 5..12, 14, 16, 18, 20
              from results/n6_mgon_P6_final.json, and the disk 0.1343093864 via a 1024-gon;
         P_7: Valtr's 2^n(3n-3)!/((n-1)!^3(2n)!) = 0.025185185 (triangle) and
              (C(2n-2,n-1)/n!)^2 = 0.033611111 (parallelogram);
         affine invariance on a NON-regular body (trapezoid vs a shear+rotation image): identical;
         4 cross-checks against convex_position.py's structurally independent monotone-chain hull
         tester (pentagon/octagon n=6, triangle/square n=7): |z| <= 2.0.
13:05  Precision anchor check at 1e10 samples (results/n6_mc_precision_1e10.json), i.e. se ~ 3e-6:
         triangle  P_6 = 0.101108954 +- 3.0e-6  (91/900,  z = -0.72)   P_7 z = -1.26
         square    P_6 = 0.122501604 +- 3.3e-6  (49/400,  z = +0.49)   P_7 z = +0.19
         disk(1024-gon), two seeds, 2e10 total: P_6 = 0.134314167 +- 2.4e-6 (exact z = +1.98)
       => the MC engine carries no detectable bias down to ~3e-6, four orders of magnitude below
       the conjecture window (0.0332).  It also sharpens the campaign's check of Marckert's
       n = 7 disk value: exact 0.0390905623 vs 0.039092252 +- 1.4e-6 here (z = +1.23; the
       previous check in LITERATURE.md was at 2e8 samples, z = +1.55).
13:06  Body catalogue `src/n6_bodies.py`: 568 area-normalised convex bodies in 10 families --
       regular m-gons m = 3..20, 24, 32, 48, 64, 128, 256, 1024 (controls); 86 vertex-pushed /
       -pulled m-gons; 62 quadrilaterals (trapezoids, skew trapezoids, kites, random); 40 random
       pentagons/hexagons; 55 Minkowski interpolations (triangle<->disk, square<->disk,
       triangle<->square, pentagon<->disk, triangle<->hexagon); 31 near-triangle bodies (1, 2 or
       3 corners truncated by t = 0.005..0.45, and triangles rounded by a disk of radius
       0.002..0.2); 62 near-disk bodies (support functions h = 1 + eps cos(k theta), k = 3..12,
       eps up to the convexity limit 1/k^2, plus rounded polygons and ellipse controls); 44
       classical smooth bodies (circular segments, sectors, stadiums, lenses, Reuleaux polygons
       of order 3, 5, 7, 9 -- verified to have constant width -- and cone/"ice-cream" shapes);
       162 random convex polygons with 4..12 vertices from three point processes.
       Smooth bodies use 256-gon (disk: 1024-gon) approximations; the polygonal bias is
       a_4/m^4 ~ 3.27/256^4 = 7.6e-10, far below the MC error.
14:35  Main scan DONE (src/n6_scan.py, 568 bodies x 4e8 samples, 88 min,
       results/n6_extremal_scan.json).  **NO VIOLATIONS.**  Over all 568 bodies the extreme
       z-scores are  min (P_6 - 91/900)/se = -1.26  (attained by a body that IS a triangle) and
       min (P_6(disk) - P_6)/se = -1.17  (the regular 32-gon, an inner approximation to the
       disk); at n = 7, -2.22 and -2.03.  Threshold for a candidate counterexample was -3.
       Control: the 27 catalogue bodies that happen to be affinely triangles scatter around
       91/900 with chi^2/df = 0.70 and mean z = +0.01.
       Regular m-gons reproduce the exact monotone table to within MC error out to m = 17; the
       five apparent decreases beyond that are all sub-sigma (the exact spacing 3.27/m^4 falls
       below the 1.7e-5 sampling error there).  The Minkowski families ending at the DISK
       (tri/sq/pent -> disk) increase strictly at all 30 steps; the two ending at a POLYGON
       (tri->sq, tri->hex) peak ABOVE both endpoints (0.128271 at t=0.7, 0.132864 at t=0.8) --
       correct, since a Minkowski sum of a triangle and a square is a heptagon.
       corr(P_6, P_7) = 0.99919 over the 568 bodies (rank corr 0.99974).
16:00  LOCAL ANALYSIS (results/n6_local_perturbation.json) -- the sharpest evidence, since a scan
       can only fail to find a counterexample.  Common random numbers give a 14x variance
       reduction on differences, so shape derivatives are measurable directly.
       DISK: fitting dP_6 = -c2 eps^2 - c3 eps^3 along h = 1 + eps cos(k theta),
         c2(k=2) = +0.0005 +- 0.0004  (ZERO: the affine/ellipse null direction)
         c2(k=3) = -0.3765 +- 0.0016   (235 sigma)      c2(k=4) = -0.5476 +- 0.0075   (73 s)
         c2(k=5) = -0.6608 +- 0.0145   (46 s)           c2(k=6) = -0.7151 +- 0.0268   (27 s)
         c2(k=8) = -0.6759 +- 0.0631   (11 s)
       => the disk is a STRICT LOCAL MAXIMISER, second variation negative definite off the
       affine direction.  Self-check that the code could not fake: the k=2 response is QUARTIC,
       dP_6 = -7.4e-7, -1.03e-5, -8.2e-5 at eps = 0.0625, 0.125, 0.225 (ratios 13.9, 7.9 vs
       2^4 = 16, 1.8^4 = 10.5) -- exactly as required, since 1 + eps cos 2theta is an ellipse to
       first order but not to second.
       TRIANGLE: cutting kc = 1,2,3 corners back a fraction t (removing area kc t^2) raises P_6
       in 14 of 15 cases (10 at >= 3 sigma), with dP_6 ~ 0.30 x (removed area fraction) and the
       SAME coefficient for 1, 2, 3 corners -- corners contribute additively.  So the triangle is
       a STRICT LOCAL MINIMISER and is a CORNER of shape space: the one-sided derivative is
       strictly positive, not zero.
       Recorded subtlety: rounding instead (T (+) rB) gives dP_6 ~ r^1.65 (exponent 1.5-1.9),
       i.e. VANISHING first derivative.  Explanation: the inner parallel bodies T (-) rB of a
       triangle are again triangles, so r -> P_6(T (+) rB) is constant 91/900 for r <= 0 and
       >= 91/900 for r > 0, forcing a zero one-sided derivative if C^1.  Exact exponent
       UNVERIFIED.
17:20  NELDER-MEAD SEARCHES (results/n6_search_kgon.json, n6_search_support.json), affine group
       quotiented out by fixing three vertices / dropping the k=1 support mode.
       MIN over k-gons, k = 4..8, two starts each: ALL TEN runs collapse to the triangle,
       P_6 = 0.1011079 +- 9.5e-6 = 91/900 to 3e-6.
       MAX over k-gons: every run converges to the REGULAR k-gon --
         k=4 -> 0.1224663/0.1224737 (49/400 = 0.1225), k=5 -> 0.1292348/0.1292291 (0.1292484),
         k=6 -> 0.1318433/0.1317837 (0.1318298), k=7 -> 0.1329543/0.1329135 (0.1329615),
         k=8 -> 0.1335092/0.1333223 (0.1335163).
       MAX over support-function Fourier coefficients (K = 4, 6, 8): the highest value in any of
       ~6600 evaluations was 0.1342343, i.e. 7.5e-5 BELOW the disk.
17:55  THE WINNER'S-CURSE CONTROL (src/n6_disk_climb.py + n6_climb_recheck.py).  Nelder-Mead
       started EXACTLY at the disk, K = 3, 5, 8, ~450 shapes each at 4e7 under CRN, reported
       best-of-run gains over the disk's same-seed value of +3.64e-5, +4.26e-5, +3.90e-5
       (2.5-3 sigma of the CRN difference error) -- i.e. three apparent counterexamples.
       Re-measured with FRESH seeds at 1.2e10 samples: 0.134306034, 0.134308587, 0.134305635
       +- 3.1e-6, i.e. -1.08, -0.26, -1.21 sigma BELOW the exact disk value.  ALL THREE
       EVAPORATED.  Exactly the outcome the protocol exists to produce.
18:35  EXACT-MACHINERY CONFIRMATION on NON-REGULAR bodies (src/n6_exact_confirm.py,
       results/n6_exact_confirm.json).  The first-principles route (width + two-chord) had only
       ever been checked on triangles/squares/regular m-gons.  Seven irregular bodies, exact vs
       4e9-sample MC, ALL OK, |z| <= 1.19:
         trapezoid(0.4)   0.11524543345034665   MC 0.115248497 +- 5.0e-6   z=+0.61
         rhombus          0.12249999999992331   (= 49/400 to 8e-13)        z=-0.20
         tri, 1 corner cut 0.05  0.10178878869638569                       z=-1.07
         tri, 3 corners cut 0.2  0.12096187020153892                       z=-0.21
         pentagon pushed 1.6x    0.12751177127659825                       z=+0.74
         random quadrilateral    0.10298445217813361                       z=+0.51
         random pentagon (quad)  0.11049377927283081                       z=-1.19
       Two structurally disjoint pipelines (a 4-D line-space integral in 80-bit arithmetic vs a
       hull test on 4e9 samples) agreeing to 5e-6 on seven irregular shapes -- and seven new
       exact P_6 values for non-regular bodies as a by-product.
19:40  PART 1 WRITE-UP: docs/N6_EXTREMAL.md.  VERDICT: the n >= 6 extremal conjecture is
       CONSISTENT with everything measured, and the two conjectured extremisers are confirmed to
       be STRICT LOCAL extremisers with measured shape derivatives.  No counterexample; no claim
       beyond consistency.

## 2026-08-19 (server) — TASK_N6_EXTREMAL Part 2: the R^d analogues

19:55  PART 2.1 SOLVED, and better than asked.  The task asked for the Dehn-Sommerville analogue
       of E[A_4] = 2E[A_3] (R^2) / E[A_5] = (5/2)E[A_4] (R^3) and for the dimensions in which the
       "one relation collapses the vertex-count moment" phenomenon survives (the R^3 note
       conjectured ODD dimensions).  Both halves have clean answers, and they point OPPOSITE ways:

       (a) THE IDENTITY IS DIMENSION-FREE, and is POINTWISE, not merely in expectation.  For any
       d+2 points in general position in R^d,
           vol conv(P_1..P_{d+2}) = (1/2) sum_{i=1}^{d+2} vol conv(P_1,..,^P_i,..,P_{d+2}).
       Proof: d+2 points in general position form a CIRCUIT -- the affine dependence
       sum lambda_i P_i = 0, sum lambda_i = 0 is unique up to scale with all lambda_i nonzero.
       Let (P,N) be the Radon partition (signs of lambda).  A circuit has exactly two
       triangulations, {conv(all minus i) : i in P} and {conv(all minus i) : i in N}; each covers
       the hull with disjoint interiors, so sum_{i in P} D_i = sum_{i in N} D_i = vol conv(all),
       and adding gives sum over ALL i = 2 vol conv(all).  Taking expectations for i.i.d.
       absolutely continuous points:
           **E[A_{d+2}] = ((d+2)/2) E[A_{d+1}]   in R^d, for EVERY d >= 1.**
       d=2 -> E[A_4] = 2E[A_3]; d=3 -> E[A_5] = (5/2)E[A_4]; both recovered.  No Efron, no
       Renyi-Sulanke, no Euler needed.  (Elementary; the two-triangulations-of-a-circuit fact is
       standard, so treat the identity as folklore-or-known, not new.)

       (b) THE DEHN-SOMMERVILLE ROUTE, by contrast, DIES AT d = 4 -- in even AND odd dimensions.
       The route needs f_{d-1} = alpha f_0 + beta for simplicial d-polytopes (because
       E[f_{d-1}(n)] = C(n,d) E[(1-c)^{n-d} + c^{n-d}] is a c-moment while E[f_0(n)] = n(1 -
       E[A_{n-1}]) is Efron).  The f-vectors of simplicial d-polytopes have floor(d/2) degrees of
       freedom modulo Dehn-Sommerville, which equals 1 exactly for d = 2, 3.  For d >= 4, f_0 and
       f_{d-1} are Dehn-Sommerville-INDEPENDENT (in d=4 the relations are f_2 = 2f_3 and
       f_1 = f_0 + f_3, and f_1 is not a c-moment).  So the R^3 note's "odd dimensions" guess is
       wrong both ways: the derivation stops at d=3, the identity never stops.
       Also: eliminating m_2 = E[c^2] between n = d+2 and n = d+3 in the general
       (alpha, beta) route reproduces exactly E[A_{d+2}] = ((d+2)/2)E[A_{d+1}] for d = 2 and 3 --
       the constants conspire, as they must.

19:55  VERIFICATION (src/rd_identity.py, results/rd_identity.json):
       * POINTWISE identity, 400 random gaussian configurations per dimension, hull volume by
         Qhull vs (1/2) sum of leave-one-out simplex determinants: worst relative error
         9.0e-16 (d=2), 1.3e-15 (d=3), 7.1e-15 (d=4), 2.5e-15 (d=5), 3.5e-15 (d=6), 4.3e-15
         (d=7).  Both Radon triangulations separately reproduce the hull volume to <= 7.6e-15.
       * DEHN-SOMMERVILLE independence, 3000 random polytopes per dimension: f_{d-1} IS a
         function of f_0 for d = 2, 3 and is NOT for d = 4, 5, 6 (e.g. in R^4 f_0 = 6 occurs with
         f_3 in {8, 9}; in R^5 f_0 = 7 with f_4 in {10, 12}; in R^6 f_0 = 8 with f_5 in
         {12, 15, 16}) -- exactly the floor(d/2) = 1 criterion.
       * IN EXPECTATION, by MC with the two sides computed by INDEPENDENT code (numpy simplex
         determinants vs Qhull hull volumes), 1e5 configurations per body:
         d = 3, 4, 5 x {ball, cube, simplex, cross-polytope, cylinder, cone} = 18 checks,
         ratio E[A_{d+2}]/E[A_{d+1}] vs (d+2)/2:  ALL |z| <= 1.19.
21:00  PART 2.2: the 3-D Sylvester five-point scan (src/rd_scan.py, results/rd_scan.json).
       P_5(K) = 1 - 5 E[A_4] (five points in R^3 are in convex position iff none is in the
       tetrahedron of the other four, and at most one can), so only a 4-point determinant
       expectation is needed -- no hull code.  240 bodies at 4e7 tetrahedra (se ~ 1e-5), 55 min:
       simplices/cubes/ball/ellipsoid anchors (7 exact checks, |z| <= 1.10), Platonic solids,
       48 prisms/antiprisms, 24 bipyramids, 32 pyramids, 31 smooth (cylinders, cones, bicones,
       capsules, ball caps), 17 near-simplex, 30 hulls of 6..256 random sphere points, 48 random
       polytopes.  **NO VIOLATIONS**: min (P_5 - simplex)/se = -1.96 (a body that IS a simplex),
       min (ball - P_5)/se = -0.85 (the ball itself).
       Ordering: simplex 0.913014 < square pyramid 0.921103 < triangular prism 0.923210 <
       cone 0.923823 < bipyramid3 0.924575 < cube 0.930781 < octahedron 0.931816 <
       hexagonal prism 0.932866 < half-ball 0.932911 < cylinder 0.933309 < bicone 0.934168 <
       dodecahedron 0.935878 < icosahedron 0.936192 < capsule(L=1) 0.936682 <
       hull of 256 sphere points 0.937046 < ball 0.937071.  Affine invariance appears unbidden:
       cone-h = 0.9238478/0.9238229/0.9238711 for h = 0.25/1/4, cylinder-h = 0.9332906/0.9333088/
       0.9333046 for h = 0.25/1/8 -- constant within error, as affine images must be.
21:37  3-D LOCAL ANALYSIS (src/rd_local.py, results/rd_local.json, 4e8 tetrahedra/body, se ~ 3e-6,
       compared against the EXACT anchors so no base run is needed):
       SIMPLEX: truncating 1 or 4 corners by t = 0.05..0.35 raises P_5 in all 7 cases, at
       +25 to +3585 sigma, with dP_5 LINEAR in the removed volume fraction and the same
       coefficient (~0.26 at the smallest cut, rising as t -> 0) for 1 and 4 corners.
       => the simplex is a STRICT LOCAL MINIMISER and a corner of shape space.
       BALL: two independent one-sided perturbations both strictly decrease P_5 --
         capsule (ball (+) segment of length L): -1.14e-5, -3.12e-5, -8.09e-5, -2.10e-4 at
           L = 0.15, 0.25, 0.40, 0.70  (-3.8 to -69 sigma); dP_5/L^2 = -0.00050 at the first
           three L, i.e. QUADRATIC in L hence SUBLINEAR in the added volume -- vanishing first
           derivative, the 3-D echo of the 2-D triangle-rounding r^1.65;
         ball cut to a cap of height 1.95/1.90/1.80/1.60: -3.97e-6, -4.20e-5, -2.38e-4,
           -1.02e-3 (-1.3 to -329 sigma); here the response IS first-order in the removed volume
           (~ -0.008 per unit), which a one-sided direction at a maximum permits.
       AFFINE NULL DIRECTION confirmed to 1e-6: ellipsoids 1.05:1:1 and 2:1:0.5 give the ball's
       exact 134/143 at -0.2 and +1.0 sigma.
       => the ball is a STRICT LOCAL MAXIMISER.
21:50  PART 2 WRITE-UP: docs/RD_ANALOGUES.md.
21:45  FINAL STATE: no campaign processes running.  Deliverables of TASK_N6_EXTREMAL:
         src/mcp6.py            (--validate: 24 checks, 17 exact anchors, ALL OK)
         src/n6_bodies.py       (568-body catalogue), src/n6_scan.py
         src/n6_search.py, src/n6_disk_climb.py, src/n6_climb_recheck.py,
         src/n6_exact_confirm.py
         src/rd_identity.py, src/rd_scan.py (--resume), src/rd_local.py
         results/n6_mc_validation.json, n6_mc_precision_1e10.json, n6_disk_P7.json,
         results/n6_extremal_scan.json (568 bodies, P_6 + P_7 + gaps + z),
         results/n6_local_perturbation.json, n6_search_kgon.json, n6_search_support.json,
         results/n6_disk_climb.json, n6_climb_recheck.json, n6_exact_confirm.json,
         results/rd_identity.json, rd_scan.json (240 bodies), rd_local.json
         docs/N6_EXTREMAL.md (new), docs/RD_ANALOGUES.md (new),
         VERDICT.md addendum, README.md (contents + 3-D simplex anchor).
       CORRECTION applied before finishing: an intermediate draft treated P_7(disk) as having no
       closed form.  It does -- Marckert 2017, Table (7), 0.0390905623, already recorded in
       LITERATURE.md.  The scan's P_7 upper anchor was switched to Marckert's exact value
       (min z below the disk at n=7 moves from -1.86 to -2.03, still no violation) and the
       2e10-sample run is now reported as a 10x-sharper CONFIRMATION of it (z = +1.23), not as a
       new anchor.
