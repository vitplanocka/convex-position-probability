# TASK (server): the polygon two-chord integral -> first-principles P_6

*Written 2026-08-19 by the local session. Work autonomously in
`~/math/convex-position-probability/`. Read `docs/N6_INGREDIENTS.md` FIRST (it has the full
derivation, the verified ingredient table, and the two-chord method already validated on the
disk). Then do this. Same verification discipline and rules as `TASK.md`.*

## Goal

Implement the two-chord Blaschke-Petkantschin integral for `E[A_4^2 & convex]` for an arbitrary
convex POLYGON, so that `P_6(K)` is assembled from first principles (no reliance on Valtr or
Marckert). Then produce closed-form / high-precision `P_6` for regular m-gons.

## The integral (derived + disk-validated; see docs/N6_INGREDIENTS.md)

    E[A_4^2 & convex] = (3/4)(1/V^6) int_{phi1,p1} int_{phi2,p2} sin^2(phi1-phi2)
                          G(a1,b1,sigma1) G(a2,b2,sigma2) 1{X in K} dp1 dphi1 dp2 dphi2

* line i = { x . n(phi_i) = p_i }, n(phi) = (cos phi, sin phi); phi in [0,pi), p in the polygon's
  projection range on n.
* [a_i, b_i] = the chord K ^ line_i as coordinates along the line (tangent t(phi)=(-sin,cos));
  i.e. clip the polygon by the line and take the two intersection points' tangential coords.
* X = line1 ^ line2; sigma_i = X's tangential coordinate along line i. Require X in K (equivalently
  a_i <= sigma_i <= b_i for both i -- and note G vanishes at the chord ends so the boundary is
  continuous).
* G(a,b,sigma) = ((b-a)^5 - (b-sigma)^5 - (sigma-a)^5)/10.
* CRITICAL normalization: the prefactor is 1/V^6, NOT 1/V^4 (D is a physical determinant; the
  extra 1/V^2 comes from E[(D/V)^2 ...]). On the disk the wrong power shows up as an exact
  factor pi^2. See src/n6_twochord_disk.py for the working disk reference (3-D, rel. 1.5e-12).

## Steps

1. Validate the machinery on the disk first by reproducing src/n6_twochord_disk.py's result
   (2400 pi^2 + 19019)/(19200 pi^4) = 0.0228343519..., then move to polygons.
2. Implement the 4-D integral for a convex polygon. Use polygon-line clipping for [a_i,b_i] and
   for X. Reduce dimension with the body's dihedral symmetry where possible (square: restrict
   phi1 to [0,pi/2) or [0,pi/4) with the right multiplicity; regular m-gon: [0, pi/m)). Integrate
   the two angular and two offset variables; for fixed (phi1,phi2) the p1,p2 dependence is smooth
   (G is piecewise polynomial), so an adaptive or fine fixed grid works. numba/vectorized numpy
   is fine for a high-accuracy numeric value; push to >= 8 digits.
3. VALIDATE: the polygon integral must reproduce, INDEPENDENTLY of Valtr,
       E[A_4^2 & convex](triangle) = 119/4500 = 0.0264444...
       E[A_4^2 & convex](square)   = 1307/54000 = 0.0242037...
   (these are in docs/N6_INGREDIENTS.md and results/n6_ingredients.json). Only once both match to
   >= 6 digits is the route trusted.
4. Then compute E[A_4^2 & convex] for regular m-gons m = 3..12 to high precision, combine with the
   exact E[A_3^3] (src/n6_bp_polygon.py: EA3k_polygon) and E[A_5] (EA_hull_polygon) to get
   E[A_4^2] and P_6(regular m-gon). PSLQ each P_6 (m=3,4,6 should be rational; m=5,10,12 in
   Q[sqrt5]/Q[sqrt3]; use the right basis). Cross-check: m=3 -> 91/900, m=4 -> 49/400,
   m -> infinity -> the disk 0.134309386.
5. If a clean closed form P_6(regular m-gon) emerges (as it did for P_5), record it. Note that
   the extremal ordering (triangle minimises, disk maximises) should show P_6 increasing in m --
   this is a data point on the OPEN n>=6 extremal conjecture (Marckert-Rahmani proved only n=5).

## Deliverables
- `src/n6_twochord_polygon.py` (the implementation), `results/n6_twochord_polygon.json`.
- Append results + any closed form to `docs/N6_INGREDIENTS.md` and log milestones in
  `PROGRESS.md` (timestamps via `date '+%F %T'`).
- If P_6(regular m-gon) closed form found: add to the table and to `docs/N5_LANDSCAPE.md`'s sibling.
- Update `VERDICT.md` addendum: E[A_4^2 & convex] now computed independently for polygons ->
  P_6 first-principles for regular m-gons.

## Rules (unchanged)
Stay in this folder; `nice -n 10`; <= 16 threads; no git; do not touch `~/math/nonmonic-01` or
`~/math/nonmonic-cubic-lean`. Reproduce anchors through every new code path before trusting a new
number; two independent computations per headline; PSLQ only unprompted from many digits. When
done (or blocked), write it up and stop.
