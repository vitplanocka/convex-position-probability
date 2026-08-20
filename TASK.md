# TASK: Convex-position probability campaign — overnight session 1

*Briefing for the autonomous Claude Code session on the server, written 2026-08-18 23:36
by the local session. Work autonomously in this folder (`~/math/convex-position-probability/`).
Read `README.md`, `PROGRESS.md`, `RESUME.md` first, then this file. Rules at the bottom are
non-negotiable.*

## Context in one paragraph

P_K(n) = P(n i.i.d. uniform points in a planar convex body K are in convex position).
Known: Valtr's exact formulas for parallelogram/triangle, Sylvester n=4 values, Marckert
2017 recursion for the disk (no closed form claimed), Barany 1999 asymptotics, the open
extremal conjecture P_triangle(n) <= P_K(n) <= P_ellipse(n). Today (day one) the local
session built a validated pipeline (two independent MC testers + a hull-area-moment route,
all reproducing 9 anchors) and derived, symbolically, three things that need to be
nailed down, verified independently, and checked for novelty:

  (I)   E[A_4] = 2 E[A_3]  for EVERY planar convex body
        (A_k = hull area of k random points as a fraction of |K|).
  (II)  P_5(K) = 1 - 10 E[A_3(1 - A_3)] = (5/2) P_4(K) - 3/2 + 15 det(Sigma_K)/|K|^2
        for EVERY planar convex body (Sigma_K = covariance matrix of the uniform law on K).
  (III) P_5(disk) = 1 - 305/(48 pi^2) = 0.35618831227...  (special case of II).

Derivation chain (each link is classical; the composition may or may not be):
Renyi-Sulanke edge count E[N_n] = C(n,2) E[(1-c)^{n-2} + c^{n-2}] with c the area fraction
on the left of the directed chord through two random points; E[c] = 1/2 by side symmetry
gives E[N_4] = 12E[c^2], E[N_5] = 30E[c^2] - 5; Efron E[A_{n-1}] = 1 - E[N_n]/n gives (I);
the Buchta-type identity P_n = sum_j (-1)^j C(n,j) E[A_{n-j}^j] gives
P_5 = 1 - 5E[A_4] + 10E[A_3^2]; and E[A_3^2] = (3/2) det Sigma/|K|^2 because
E[det^2] = 6 det Sigma is a polynomial moment. Checks so far: (II) reproduces Valtr's
11/36 (triangle, E[A_3^2] = 1/72) and 49/144 (square, 1/96); (III) matched by direct MC
(2e8, z = -0.4) and route M (z = -0.09); sympy exact in `src/disk_n5_exact.py`.

## Mission, in priority order

### A. Make (I)-(III) bulletproof  (first 1-2 hours)
1. Run the anchors through every code path first: `python src/convex_position.py --validate`,
   `python src/route_moments.py --validate`, `python src/general_n5_identity.py`
   (numba is in `.venv`; use `../.venv/bin/python` from `src/`).
2. Write `THEOREMS.md`: full statement + proof of (I) and (II) with every step justified,
   including a clean proof of E[c] = 1/2, of E[det^2] = 6 det Sigma, and of the Buchta-type
   identity (derive it: P(N_n = m) = C(n,m) E[1_conv(m) A_m^{n-m}] then invert). Check every
   symbolic step with sympy where possible. Note explicitly what fails at n = 6
   ((1-c)^4 + c^4 does not collapse; P_6 needs E[A_4^2], E[A_3^3]).
3. Independent numerical confirmation of (I) and (II) on at least six bodies: square,
   triangle, disk, regular pentagon/hexagon/octagon, a long ellipse (trivial by affine
   invariance, but run it), a half-disk, a stadium. For each: E[A_3], E[A_4], E[A_3^2] by MC
   (>= 1e8 triples), P_4 and P_5 by the DIRECT tester (>= 1e8), and compare with (II).
   Record z-scores. Anything beyond 3 sigma is a bug until proven otherwise.
4. Novelty check (WebSearch/WebFetch are available to you). Highest priority sources:
   Buchta 2005 'An identity relating moments of functionals of convex hulls' (DCG 33:125-142)
   -- read what its applications section actually states; Marckert 2017 arXiv:1402.3512
   (does it give p_5 = 1 - 305/(48 pi^2)?); Efron 1965 Biometrika; Renyi-Sulanke 1963/64;
   Hilhorst-Calka-Schehr 2008; classical books (Kendall-Moran, Santalo, Solomon, Mathai).
   Search literal strings "305/48", "305/(48", "0.35618", "175/24", "35/(24", "twice"
   "expected area" four points three points. If `docs/lit_*.md` files exist (the local
   session's literature agents drop them there), read them first and cite them; if not,
   do your own sweep and write `LITERATURE.md` (facts tagged [verified]/[secondary]/
   [unverified], every claim with a URL). Verdict must be explicit: known / partially known
   (which part) / not found (list what was searched).

### B. Exact P_5 for regular polygons and the n = 5 extremal landscape  (main overnight job)
By (II), P_5(K) needs only P_4(K) (Sylvester's four-point problem, i.e. E[A_3]) and
det Sigma_K/|K|^2 (exact for any polygon: see `polygon_cov_det` in
`src/general_n5_identity.py`). So:
1. Compute E[A_3] EXACTLY for regular m-gons, m = 3..12 (and if feasible general m):
   symbolic integration (sympy) of the random-triangle mean area in a polygon, or reproduce
   the classical regular-polygon Sylvester values from the literature (there are published
   formulas/tables for regular polygons -- find and CITE them, then reproduce them
   independently; do not just copy). Anchor: square 11/144, triangle 1/12, hexagon and
   pentagon must match MC to >= 5 digits.
2. Tabulate exact P_4, E[A_3^2], P_5 for m = 3..12 and the disk limit; PSLQ/identify any
   surprising constants only when they come unprompted from >= 20 digits.
3. Extremal conjecture at n = 5: since P_5 = 1 - 10 E[A_3(1-A_3)], the n=5 conjecture is
   equivalent to: the functional F(K) = E[A_3(1-A_3)] is maximised by triangles and
   minimised by ellipses. Blaschke gives E[A_3] extremal in that direction, but
   E[A_3^2] = (3/2) det Sigma/|K|^2 is extremal in the OPPOSITE direction (ellipse minimises,
   triangle maximises the planar isotropic constant), so the conjecture is NOT immediate.
   (a) Test it numerically to high precision across families: regular m-gons, m-gons with
   one vertex pushed, ellipse-triangle interpolations (Minkowski combinations), stadiums,
   half-disks, random convex polygons (>= 200 random bodies), affine-normalised. Report the
   closest approaches to the bounds. (b) Look for a proof route: e.g. is F monotone under
   Steiner symmetrisation? Under 'shaking'? Is there a known inequality between E[A_3^2] and
   E[A_3]^2 that closes the gap? Try; report honestly whether anything works.
4. Also compute the n=5 landscape via the DIRECT tester for a subset, as an independent
   check of the whole (II)-based table.

### C. n = 6, disk  (only if A and B are done and written up)
P_6 = 1 - 6E[A_5] + 15E[A_4^2] - 20E[A_3^3]. E[A_5] = 1 - E[N_6]/6 and E[N_6] = 30 E[c^4]
(show this; E[c^3] is determined by E[c^2] by symmetry). For the disk, E[c^4] is a 1-D
integral (`src/disk_semi_analytic.py` has the setup). E[A_3^3] for the disk: find or
derive the exact third moment of a random triangle's area (Crofton-style; check the
literature). E[A_4^2] = E[A_4^2 1_conv] + 4E[A_3^3]: the convex-quadrilateral second moment
is the hard part -- explore, but do not sink the night into it. Compare everything with the
direct MC value 0.13432961 +- 2.4e-5 and route M 0.1343253 +- 3.9e-5, and with whatever
Marckert 2017 gives for p_6.

## Deliverables (write to this folder as you go)
- `PROGRESS.md`: append timestamped entries (`date '+%F %T'`) at every milestone.
- `THEOREMS.md`: statements + proofs + sympy checks for (I), (II), (III), and E[N_6] = 30E[c^4].
- `LITERATURE.md` (or an appended server section if the local `docs/lit_*.md` exist):
  novelty verdicts with URLs.
- `results/`: JSON per run (body, n, samples, seed, p_hat, std_err, method).
- `results/n5_landscape.json` + `docs/N5_LANDSCAPE.md`: the exact regular-polygon table and
  the numerical extremal test.
- `VERDICT.md` at the end: what is established (with the evidence table), what is new vs
  known, what failed, and a recommended next mission.

## Verification discipline (non-negotiable, inherited from the previous campaign)
- Reproduce the anchors through every new code path before believing any new number.
- Every headline number needs two structurally independent computations that agree.
- A quadrature routine's own error estimate is not evidence; agreement between independent
  routes is. Treat endpoint singularities and near-degenerate configurations as suspect.
- PSLQ/identify hits count only when unprompted from many digits; record the basis tried.
- Log unverified claims as UNVERIFIED and leave them that way; do not let them harden.
- If a claim in this brief turns out wrong, say so in PROGRESS.md and fix it everywhere.

## Rules
- Stay inside `~/math/convex-position-probability/`. Do NOT touch `~/math/nonmonic-01`
  (a live Claude session), `~/math/nonmonic-cubic-lean` (a live Lean build), or `/var/www`.
- `nice -n 10` all heavy computation; keep RAM < 8 GB; the machine has 28 cores but is
  shared -- use at most ~16 threads (`NUMBA_NUM_THREADS=16`).
- Python: `../.venv/bin/python` (numpy 2.5, scipy 1.18, sympy 1.14, mpmath 1.3, numba 0.67).
- No git operations; the operator syncs back manually.
- When A, B are done (and C if reached), write VERDICT.md, then STOP. Do not start new
  research directions. If blocked, write what is blocking in PROGRESS.md and stop.
