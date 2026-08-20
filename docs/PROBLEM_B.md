> ## !! CORRECTION ADDED 2026-08-19 00:33 BY THE OVERNIGHT SERVER SESSION !!
>
> **The problem posed below is already solved. Do not spend time on it as an open problem.**
>
> * **The n = 5 extremal conjecture is a theorem.** Marckert & Rahmani, *Around Sylvester's
>   question in the plane*, Mathematika **67** (2021) 860–884, doi:10.1112/mtk.12110, proved
>   `Q_T^5 <= Q_H^5 <= Q_D^5` for every planar convex body. So items (a) and (b) of "The
>   problem" below are both **known to be true**. Only `n >= 6` is open. Independent
>   confirmation: Marckert & Morin, [arXiv:2411.08456](https://arxiv.org/pdf/2411.08456),
>   p. 4, states verbatim: *"The case n = 5 was proved by Marckert and Rahmani [18], but
>   n >= 6 is still a conjecture (but one finds in [18] an explicit formula to compute
>   P_K(2)(n))."*  Section 5 below ("To our knowledge n = 5 is open") is therefore wrong.
> * **The "new identity" below is NOT in Marckert & Rahmani** (checked 2026-08-19: the paper
>   was obtained via the otherwise unadvertised preprint
>   [arXiv:1511.03658](https://arxiv.org/abs/1511.03658) and read in full). Their "new formula
>   for `Q_H^n`" is a different construction — condition on the sorted abscissas, reduce to the
>   probability that uniform points on `n` vertical *segments* are in convex position, and
>   evaluate that in closed form. Across the whole paper the triangle area appears only inside
>   the classical `n=4` statement, and Efron, Rényi, Sulanke, "second moment", "covariance" and
>   `E[A_3^2]` never occur. So the identity below is **not found** in the literature — which is
>   still not the same as **new**: it is Buchta's classical identity combined with
>   `E[A_4] = 2E[A_3]`, and that in turn is three classical facts away. Residual gap: only the
>   *submitted* version is open-access; the published 2021 text is paywalled.
> * **`P_5(disk) = 1 - 305/(48 pi^2)` is in the literature**: Marckert 2017,
>   [arXiv:1402.3512](https://arxiv.org/pdf/1402.3512), Table (7). Same table gives exact
>   `P_disk(6), P_disk(7), P_disk(8)`.
> * **Corrections to the numbers below.** The Monte-Carlo values quoted for the regular
>   pentagon and hexagon are now exact:
>   `P_5(pentagon) = 7/12 - 47*sqrt(5)/450 = 0.349788455683355`,
>   `P_5(hexagon)  = 1373/3888 = 0.353137860082305`, and in general
>   `P_5(regular m-gon) = 1 - 5(15 cos^2 w + 92 cos w + 76)/(36 m^2 sin^2 w)`, `w = 2 pi/m`
>   (Alikoski 1939 + the identity). See `../docs/N5_LANDSCAPE.md`.
> * **One statement below is misleading as written.** "the two terms of F are extremised by
>   the SAME bodies but in OPPOSITE directions" — both `E[A_3]` and `det(Sigma)/area^2` are
>   maximised by triangles and minimised by ellipses; what is opposite is the *sign with
>   which they enter* `P_5`. The conclusion drawn (that the conjecture does not follow from
>   the two known extremal results) is correct.
> * What *is* still open and worth the effort: **n >= 6**. See `../VERDICT.md`, "Recommended
>   next mission".
>
> Everything below this line is the original text of 2026-08-18, left unchanged.

---

# Problem B: Sylvester's five-point problem and the extremal conjecture at n = 5

## Setting

Let K be a convex body in the plane (compact, convex, nonempty interior). Draw n points
X_1, ..., X_n independently and uniformly at random in K. They are "in convex position" if
every X_i is a vertex of the convex hull conv{X_1, ..., X_n}. Write

    P_n(K) = P( X_1, ..., X_n are in convex position ).

P_n(K) is invariant under affine maps of K (so "square" = every parallelogram, "disk" =
every ellipse). Write A_3 = area(conv{X_1, X_2, X_3}) / area(K), the normalised area of a
random triangle in K.

## Known facts (all classical, all safe to use)

1. Sylvester's four-point problem: P_4(K) = 1 - 4 E[A_3]  (four points fail to be in convex
   position iff one lies in the triangle of the other three).
   Values: triangle 2/3, parallelogram 25/36, disk 1 - 35/(12 pi^2) = 0.704480.
2. Blaschke (1917): among all planar convex bodies, E[A_3] is maximised by triangles and
   minimised by ellipses; equivalently 2/3 <= P_4(K) <= 1 - 35/(12 pi^2), with equality
   only for triangles / ellipses.
3. Valtr (1995, 1996): exact P_n for parallelograms, [C(2n-2,n-1)/n!]^2 (n=5: 49/144), and
   for triangles, 2^n (3n-3)! / [(2n)! ((n-1)!)^3] (n=5: 11/36).
4. Second moment of the random triangle area is a polynomial moment:
   E[A_3^2] = (3/2) det(Sigma_K) / area(K)^2, where Sigma_K is the covariance matrix of the
   uniform distribution on K. (Because E[det^2] = 6 det Sigma for three i.i.d. points with
   covariance Sigma, and the triangle area is |det|/2.) In the plane, det(Sigma_K)/area(K)^2
   is the square of the isotropic constant of K; it is minimised by ellipses (Blaschke) and
   maximised by triangles (known in the plane). Values: triangle 1/72, parallelogram 1/96,
   disk 3/(32 pi^2).
5. Open conjecture (attributed to Barany and others): for every n >= 4 and every planar
   convex body K,  P_n(triangle) <= P_n(K) <= P_n(ellipse). Blaschke's theorem is the case
   n = 4. To our knowledge n = 5 is open.

## New identity (derived 2026-08-18; verified on triangle, parallelogram, disk; please re-derive)

For every planar convex body K,

    P_5(K) = 1 - 10 E[ A_3 (1 - A_3) ]  =  (5/2) P_4(K) - 3/2 + 15 det(Sigma_K)/area(K)^2.

Derivation sketch: (i) Renyi-Sulanke edge count: E[#vertices of conv of n points]
= C(n,2) E[(1-c)^{n-2} + c^{n-2}], where c is the fraction of area(K) on the left of the
directed chord through two random points; E[c] = 1/2 by side symmetry, so E[N_4] = 12 E[c^2],
E[N_5] = 30 E[c^2] - 5. (ii) Efron: E[A_{n-1}] = 1 - E[N_n]/n (A_k = normalised hull area of
k points), hence E[A_4] = 2 E[A_3] for every K. (iii) Buchta-type identity from
P(N_n = m) = C(n,m) E[1{first m points in convex position} A_m^{n-m}]:
P_5 = 1 - 5 E[A_4] + 10 E[A_3^2]. Combine.
Checks: triangle gives (5/2)(2/3) - 3/2 + 10/72 = 11/36; parallelogram gives
(5/2)(25/36) - 3/2 + 10/96 = 49/144 (both agree with Valtr);
disk -> 1 - 305/(48 pi^2) = 0.3561883, matched by Monte Carlo (2e8 samples, z = -0.4).

## The problem

Because of the identity, the n = 5 case of the extremal conjecture is EQUIVALENT to:

    Among all planar convex bodies K, the functional
        F(K) = E[ A_3 (1 - A_3) ] = E[A_3] - (3/2) det(Sigma_K)/area(K)^2
    is maximised by triangles and minimised by ellipses.

Numerically: F(triangle) = 1/12 - 1/72 = 5/72 = 0.069444 (P_5 = 11/36 = 0.305556);
F(square) = 11/144 - 1/96 = 19/288 = 0.065972 (P_5 = 49/144 = 0.340278);
F(disk) = 35/(48 pi^2) - 3/(32 pi^2) = 61/(96 pi^2) = 0.064381 (P_5 = 0.356188).
Regular pentagon (Monte Carlo): P_5 ~ 0.3495; regular hexagon: P_5 ~ 0.3529.

The difficulty: the two terms of F are extremised by the SAME bodies but in OPPOSITE
directions -- E[A_3] is largest for triangles / smallest for ellipses (Blaschke), while
det(Sigma)/area^2 is ALSO largest for triangles / smallest for ellipses -- and they enter F
with opposite signs. So the conjecture at n = 5 does not follow from the two known extremal
results; it needs a genuinely finer inequality: roughly, that E[A_3] varies "faster" than
(3/2) det(Sigma)/area^2 across convex bodies.

Asked:
(a) Prove (or disprove) that F(K) <= F(triangle) = 5/72 for every planar convex body K, with
    equality only for triangles. [Lower bound of P_5 by the triangle.]
(b) Prove (or disprove) that F(K) >= F(ellipse) = 61/(96 pi^2) for every planar convex
    body K, with equality only for ellipses. [Upper bound of P_5 by the ellipse.]
Either half alone is of interest. Natural tools: Steiner symmetrisation and "shaking"
(Blaschke's original proof of the n = 4 case shows E[A_3] does not increase under Steiner
symmetrisation -- does F behave monotonically under symmetrisation or under Blaschke's
shaking of polygons?); the affine-invariant formulations (E[A_3] as an integral of chord
powers via Crofton / Blaschke-Petkantschin: E[A_3] and E[c^2] are both integrals over the
line space of chord-length powers times cap-area functions); known inequalities relating
random-triangle moments to the isotropic constant. Also useful: reduce to polygons by
approximation and prove a local (vertex-moving) monotonicity, or find a counterexample body
by optimisation. Please state clearly what is proved, what is checked numerically, and what
remains open.

## Sanity anchors for any computation
P_5(triangle) = 11/36; P_5(square) = 49/144; P_5(disk) = 1 - 305/(48 pi^2);
E[A_3]: triangle 1/12, square 11/144, disk 35/(48 pi^2);
E[A_3^2]: triangle 1/72, square 1/96, disk 3/(32 pi^2);
E[A_4] = 2 E[A_3] in every case.


---
## RESOLUTION (2026-08-19 06:30, integrated by the local session)

Both (a) and (b) are TRUE (already Marckert-Rahmani 2021, Theorems 1 and 3). An external model
session answering this brief supplied a short proof through the identity: a fibrewise "section
lemma" (flat plateau of the noise density on |t| <= r; exact gap formulas for Steiner
symmetrisation and shaking; fibre ratio Delta E[A^2]/Delta E[A] = b_2/V < 1). Written up and
numerically verified in [`N5_PROOF.md`](N5_PROOF.md).
