# TASK (server): attack the open n >= 6 extremal conjecture on non-regular bodies, + R^d

*Written 2026-08-19 by the local session. Work autonomously in
`~/math/convex-position-probability/`. Read `docs/N6_LANDSCAPE.md` and `docs/N6_INGREDIENTS.md`
FIRST (the n=6 machinery: width route for E[A_3^k]/E[A_5] in `src/n6_bp_polygon.py`, two-chord
route for E[A_4^2&conv] in `src/n6_twochord_polygon.py`, assembly `P_6 = 1 - 6E[A_5] +
15E[A_4^2&conv] + 40E[A_3^3]`). Same verification discipline and rules as `TASK.md`.*

## Background / status

The exact P_6 machinery now works for any convex polygon; `src/n6_p6_body.py` (local) assembles
P_6(V) for an arbitrary polygon and reproduces 91/900 for the right triangle. But the EXACT route
is ~15-60 s/body and too slow for a broad shape scan. The extremal conjecture window is WIDE
(P_6(triangle) = 91/900 = 0.1011111, P_6(disk) = 0.1343093864), so direct convex-position Monte
Carlo at 1e8-1e9 samples (4-5 digits) is the right tool for the SCAN, with the exact machinery
reserved to confirm any near-extremal body.

## Part 1 -- the n >= 6 extremal scan (main job)

CONJECTURE (open for n >= 6; Blaschke n=4, Marckert-Rahmani n=5):
    P_n(triangle) <= P_n(K) <= P_n(ellipse)  for every convex body K.
Test it at n = 6 (and n = 7) numerically to high confidence over a broad, affinely-normalised
family of NON-regular bodies. Affine invariance => all triangles give 91/900 and all
parallelograms 49/400, so a meaningful test needs bodies not affine-equivalent to a regular one.

1. Build a fast direct-MC `P_n(K)` for an arbitrary convex polygon: uniform sampling by
   area-weighted fan triangulation; convex position via the exact orientation test (a point is a
   non-vertex iff it lies in a triangle of three others -- reuse `convex_position.py`'s tester, or
   a numba hull-size==n test). Vectorise / numba; `nice -n 10`, <= 16 threads. Validate it
   reproduces 91/900, 49/400 and the exact regular-m-gon P_6 from `results/n6_mgon_P6_final.json`
   to within MC error before trusting it.
2. Scan a large family (aim >= 300 bodies), all area-normalised, e.g.:
   - regular m-gons m = 3..20 (controls; must be monotone increasing to the disk);
   - general trapezoids, kites, and other irregular quadrilaterals/pentagons/hexagons;
   - one-vertex-pushed and one-vertex-pulled regular m-gons;
   - ellipse<->triangle and polygon<->disk Minkowski interpolations (1-parameter families);
   - half-disk, circular segments, stadium (capsule), lens, Reuleaux-triangle approximation;
   - a few hundred RANDOM convex polygons (random points -> hull), 4..12 vertices.
   For each: P_6 (and P_7) by MC at >= 2e8 samples; record P_6, its std error, and the signed
   gaps to 91/900 and to the disk value.
3. Then push the SEARCH for extremes: minimise and maximise P_6 over polygon shape (Nelder-Mead
   or CMA over vertex coordinates, affine-normalised) starting from several seeds. The conjecture
   predicts min -> triangle (91/900), max -> disk (0.1343094).
   - Report the closest approaches to each bound and the shapes achieving them.
   - **Any body with P_6 < 91/900 - 3 sigma or > disk + 3 sigma is a potential COUNTEREXAMPLE**:
     re-run it at 1e10 samples AND confirm with the exact two-chord+width machinery
     (`src/n6_p6_body.py` pattern) before making any claim. A refuted "violation" is the expected
     outcome (the conjecture is widely believed); a surviving one would be a major result -- treat
     with extreme skepticism and triple-check.
4. Deliverable: `results/n6_extremal_scan.json` (all bodies, P_6, P_7, gaps), `docs/N6_EXTREMAL.md`
   (method, the min/max found, closest approaches, verdict: conjecture consistent / violated, with
   the evidence). Note whether the ORDERING among regular polygons and the tested families is
   monotone as the conjecture's proof strategy (Steiner symmetrisation) would predict.

## Part 2 -- R^d analogues

1. The identity family in R^d. In R^3 you already proved `E[A_5] = (5/2) E[A_4]` from Euler
   `V = 2 + F/2` (simplicial). Generalise: use the Dehn-Sommerville relations to get the analogue
   of `E[A_4] = 2 E[A_3]` / `E[A_5] = (5/2)E[A_4]` in R^d. Determine in which dimensions the
   "one relation collapses the vertex-count moment" phenomenon survives (the R^3 note conjectured
   ODD dimensions). State and prove what you can; check numerically in R^3 and R^4 (MC, exact
   orientation predicates via Caratheodory -- a point is a non-vertex iff it is in the hull of
   d+1 others).
2. The 3-D Sylvester 5-point extremal question (open analogue of Blaschke): is
   `P(d+2 points in convex position)` minimised by the simplex and maximised by the ball in R^d?
   For d = 3, MC-scan P(5 points in convex position) = 1 - 5 E[vol tetra]/|K| over 3-D bodies
   (regular simplex, cube, ball, regular octahedron, cylinder, cone, bipyramid, random convex
   polytopes). Anchors (exact): ball 134/143 = 0.937063, cube 0.930786 (Zinani). Report whether
   the simplex is the min and the ball the max, and the closest approaches.
3. Deliverable: `docs/RD_ANALOGUES.md`, `results/rd_scan.json`, with the identity statements
   (proved / conjectured), the R^3/R^4 numeric checks, and the 3-D 5-point extremal data.

## Rules (unchanged)
Stay in this folder; `nice -n 10`; <= 16 threads; no git; do not touch `~/math/nonmonic-01` or
`~/math/nonmonic-cubic-lean`. Reproduce anchors through every new code path before trusting a new
number; two independent computations per headline; MC "violations" must be reproduced at higher
sample count AND confirmed by the exact machinery before any claim; PSLQ only unprompted from many
digits; log unverified claims as UNVERIFIED. Do Part 1 fully (write it up) before Part 2. When done
write `VERDICT` updates and stop.
