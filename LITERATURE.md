# LITERATURE — verified facts, sources, and novelty verdicts

*Compiled by the overnight server session, 2026-08-19. The local session's `docs/lit_*.md`
were never synced to the server (the `docs/` folder is empty here), so this is an
independent sweep. Tags: **[verified]** = we reproduced the claim ourselves with an
independent computation; **[secondary]** = read in a reliable secondary source (MathWorld,
a citing paper) but the primary source was not opened; **[unread]** = we could not obtain
the source. Every claim carries a URL.*

## 1. Classical exact values (all [verified] by our own code)

| fact | source | our check |
|---|---|---|
| `P_triangle(4) = 2/3`, `P_parallelogram(4) = 25/36` | Sylvester's four-point problem; [MathWorld](https://mathworld.wolfram.com/SylvestersFour-PointProblem.html), [Wikipedia](https://en.wikipedia.org/wiki/Sylvester%27s_four_point_problem) | exact rational route P; MC |
| `P_disk(4) = 1 - 35/(12 pi^2) = 0.70447988...` | Woolhouse/Crofton; Blaschke 1917 extremality; [MathWorld](https://mathworld.wolfram.com/SylvestersFour-PointProblem.html) | route P Richardson to 2e-16; MC |
| `P_parallelogram(n) = [C(2n-2,n-1)/n!]^2` | Valtr, *Probability that n random points are in convex position*, Discrete Comput. Geom. 13 (1995) 637–643 | `src/anchors.py`; MC n=4..10 |
| `P_triangle(n) = 2^n (3n-3)! / [(2n)! ((n-1)!)^3]` | Valtr, Combinatorica 16 (1996) 567–573 | `src/anchors.py`; MC n=4..10 |
| `P_cube3D(5) = 1 - 5(3977/216000 - pi^2/2160)`, `P_ball3D(5) = 134/143` | Zinani 2003; classical `9/715` | MC 2e7, both testers, z<1 |

## 2. Alikoski's regular-polygon formula — the key input for section B

> "The mean triangle area of a triangle picked inside a regular n-gon of **unit area** is
> `A_n = (9 cos^2 w + 52 cos w + 44) / (36 n^2 sin^2 w)`, where `w = 2 pi/n` (Alikoski 1939)."

Source: [MathWorld, *Polygon Triangle Picking*](https://mathworld.wolfram.com/PolygonTrianglePicking.html)
**[secondary]** for the attribution — the primary is
H. A. Alikoski, *Über das Sylvestersche Vierpunktproblem*, Ann. Acad. Sci. Fenn. **51** No. 7
(1939) 1–10, which we could not obtain **[unread]**. MathWorld's table of values
`n = 3,4,5,6,8`: `1/12, 11/144, (9+2 sqrt 5)/180, 289/3888, (97+52 sqrt 2)/2304`.

**[verified], strongly.** Our route P (Proposition 6.1 of `THEOREMS.md`), which shares no
input with Alikoski, reproduces
* all five tabulated algebraic values exactly (`1/12` and `11/144` by exact `Fraction`
  arithmetic; `289/3888` likewise via an affinely-regular rational hexagon);
* the formula to **50 decimal digits** for `m = 3..24, 30, 40, 60, 100`
  (max residual `3.4e-49`, `results/B1_mgon_exact_table.json`; the built-in checks `T_0 = 2` and `E[N_3] = 3` hold to `1.9e-48`, and the `P_5` closed form of §E4 reproduces to `3.3e-51`);
* the `m -> infinity` limit `35/(48 pi^2)`, matched to `6e-17`.
Internal exactness checks `T_0 = 2` and `E[N_3] = 3` hold identically throughout.

Also **[verified]**: Buchta's exact expected hull areas
(square `11/144, 11/72, 79/360, 199/720`; triangle `1/12, 1/6, 43/180, 3/10, 197/560`;
[Buchta 1984, quoted here](https://gdr-geostoch.math.cnrs.fr/workshop_Rouen/SlidesRouen2012/Friday/Vortrag_Buchta.pdf))
— route P returns **every one of these as an exact rational**, `E[A_3]` through `E[A_7]`.

## 3. The disk for `n >= 5` — settles our Corollary 5.1

Marckert, *Probability that n random points in a disk are in convex position*, Braz. J.
Probab. Stat. 31 (2017) 320–337; preprint [arXiv:1402.3512](https://arxiv.org/pdf/1402.3512).
We read the PDF. Its Table (7) gives, **verbatim**, `1 - P^D_n` for `n = 4..8`:

    n:        4              5               6                            7                             8
    1-P^D_n:  35/(12 pi^2)   305/(48 pi^2)   (146400 pi^2 - 473473)/(11520 pi^4)
                                             (512400 pi^2 - 2900611)/(23040 pi^4)
                                             (62664108221 + 1721664000 pi^4 - 18670881600 pi^2)/(48384000 pi^6)

and also `P^D_{5,3} = 15/(16 pi^2)`, `P^D_{5,4} = 65/(12 pi^2)`, etc.

**Consequence: `P_disk(5) = 1 - 305/(48 pi^2)` is KNOWN, not new.** The local session flagged
this as "NOVELTY UNKNOWN"; the verdict is now: **known, Marckert 2017**.

**[verified]** — our Monte Carlo (2e8 samples, direct tester) against Marckert's exact values:

| n | Marckert exact | our direct MC | z |
|---|---|---|---|
| 5 | 0.3561883123 | 0.356188180 ± 3.4e-5 | −0.00 |
| 6 | 0.1343093864 | 0.134329610 ± 2.4e-5 | +0.84 |
| 7 | 0.0390905623 | 0.039112300 ± 1.4e-5 | +1.55 |
| 8 | 0.0090342622 | 0.009042290 ± 6.7e-6 | +1.20 |

This simultaneously validates our MC pipeline and confirms Marckert's table.
(It also means TASK.md section C's target for `n=6` in the disk is already known exactly.)

## 4. The extremal conjecture — n = 5 IS ALREADY A THEOREM

Marckert & Rahmani, *Around Sylvester's question in the plane*, Mathematika **67** (2021)
860–884, [doi:10.1112/mtk.12110](https://londmathsoc.onlinelibrary.wiley.com/doi/10.1112/mtk.12110);
preprint **[arXiv:1511.03658](https://arxiv.org/abs/1511.03658)**. **[verified]** — obtained and
read in full (11 075 words). Wiley (402) and HAL (anti-bot) stayed blocked and the author's own
publication page lists only the HAL link; the arXiv id came from the Semantic Scholar Graph API
on the DOI. Unpaywall and OpenAlex both confirm every OA copy (arXiv *and* HAL) is the
**submitted** version, so this is the best obtainable text and the published version stays
paywalled. They prove `Q_T^5 <= Q_H^5 <= Q_D^5`, and their "new formula for `Q_H^n`" is:

    Q^{N+2}_H = int_{ND_{N+2}} <V_H(x_0), ..., V_H(x_{N+1})> f_H(x[N+2]) dx,

i.e. condition on the sorted abscissas, so that the problem reduces to `<S_0,...,S_{N+1}>`, the
probability that independent uniform points on `N+2` given **segments** are in convex position;
the new content is a closed rational-fraction evaluation of that (Prop. 17), expanded over
triangulations (Thm 20) and permutations (Lem. 23). The point is that the density
`f_H = (N+2)! prod W_H(x_j)` depends only on the width function, hence is invariant under
Steiner symmetrization and shaking — Blaschke's n=4 method, pushed to n=5.

Their Theorem 1 as typeset reads `Q^5_D = 1 - 305/(48 pi)^2`, which evaluates to 0.98659 and is
impossible; the intended reading is `1 - 305/(48 pi^2) = 0.35618831...` (and likewise
`Q^4_D = 1 - 35/(12 pi^2)`). Bracket typos in the arXiv version — but this is a **second
independent source** for `P_5(disk)`, agreeing with Marckert 2017 Table (7).

Confirmed **[secondary]** by a direct quotation in a later paper we did read in full —
Marckert & Morin, *The Sylvester question in R^d: convex sets with a flat floor*,
[arXiv:2411.08456](https://arxiv.org/pdf/2411.08456), p. 4:

> "The case n = 5 was proved by Marckert and Rahmani [18], but n ≥ 6 is still a conjecture
> (but one finds in [18] an explicit formula to compute P_K(2)(n))."

**So TASK.md section B3 ("test the n=5 conjecture, look for a proof route") is a solved
problem.** What we did instead — a high-precision *independent verification* of the theorem
over ~500 exactly-computed convex bodies plus a direct numerical optimisation — is a check
of a published result, not a discovery. Reported as such in `docs/N5_LANDSCAPE.md`.

Related and read: Bárány, *Sylvester's question: the probability that n points are in convex
position*, Ann. Probab. 27 (1999) 2020–2034
([Project Euclid](https://projecteuclid.org/journals/annals-of-probability/volume-27/issue-4/Sylvesters-Question--The-Probability-That-n-Points-arein-Convex/10.1214/aop/1022874826.full));
Morin, regular-`kappa`-gon and general-polygon asymptotics
([arXiv:2401.16207](https://arxiv.org/pdf/2401.16207), [arXiv:2410.11706](https://arxiv.org/abs/2410.11706))
— both **asymptotic** in `n`, no small-`n` closed forms.

## 5. Novelty verdicts for this session's results

| result | verdict | basis |
|---|---|---|
| **(III)** `P_disk(5) = 1 - 305/(48 pi^2)` | **KNOWN** | Marckert 2017 Table (7), read verbatim |
| **n=5 extremal conjecture** | **KNOWN (a theorem)** | Marckert–Rahmani 2021, quoted by Marckert–Morin 2024 |
| **(I)** `E[A_4] = 2 E[A_3]` for every abs. continuous planar law | **NOT FOUND as a stated result** — but it is one step from three classical facts, is *visible* in Buchta's published tables (`11/144 -> 11/72`, `1/12 -> 1/6`), and its Gaussian instance is noted in the literature (below). Treat as folklore-or-known; **do not claim as new**. | see search log |
| **(I')** `E[A_6] = 3E[A_5] - 5E[A_3]`, **(I'')** `E[A_8] = 4E[A_7] - 14E[A_5] + 28E[A_3]`, and the general "one new parameter every other step" structure | **NOT FOUND** | see search log |
| **(II)** `P_5 = 1 - 10 E[A_3(1-A_3)]` (and `= (5/2)P_4 - 3/2 + 15 det Sigma/|K|^2`) | **NOT FOUND** (upgraded from UNRESOLVED on 2026-08-19 once Marckert–Rahmani was obtained and read). Their "new formula" is a different object entirely — see §4. In the whole paper the triangle area appears only inside the classical `n=4` statement, and there is no occurrence of Efron, Rényi, Sulanke, "second moment", "covariance", "isotropic", `E[A_3^2]`, or anything of the form `Q^5_H = 1 - 10 E[...]`. Same caveat as (I): it is Buchta's classical identity plus (I), so **do not claim as new**. | §4; `PROGRESS.md` 06:34 |
| **Corollary 5.2**, `P_5(regular m-gon) = 1 - 5(15cos^2 w + 92 cos w + 76)/(36 m^2 sin^2 w)` | **NOT FOUND** in any source we searched; but it is a one-line corollary of Alikoski (1939) + (II), so the novelty is only as strong as (II)'s. Report as "explicit closed form, apparently not tabulated anywhere; derivation elementary given Alikoski". | see search log |
| **(I_3)** `E[A_5] = (5/2) E[A_4]` in `R^3` (and `E[A_7] = (7/2)E[A_6] - (35/4)E[A_4]`), proved from the Renyi-Sulanke facet count + `V = 2 + F/2` for simplicial 3-polytopes + Efron | **NOT FOUND** | searched "expected volume convex hull five random points three dimensions 9/715 Efron identity"; read Buchta, *The duality of the volumes and the numbers of vertices of random polytopes*, [arXiv:2107.05402](https://arxiv.org/pdf/2107.05402) (about a different transformation of the Efron–Buchta right-hand side, no such relation) and Kabluchko–Last–Zaporozhets, [arXiv:1603.01357](https://arxiv.org/abs/1603.01357) (combinatorial inclusion–exclusion + Euler relation; does not mention Efron or Sylvester at all). Same caveat as (I): one step from classical ingredients. |
| **Route P** (edge-pair cell decomposition of line space giving exact `E[N_n]`, `E[A_{n-1}]` for any convex polygon) | **method presumably known** — Santaló's boundary parametrisation of `dp dtheta` is classical and Alikoski/Buchta evidently did polygon computations. Our contribution is a clean, fully verified implementation, not a new idea. | — |

**The one relevant near-hit for (I).** Finch, *Capturing, Ordering and Gaussianity in 2D*,
[arXiv:1601.04937](https://arxiv.org/pdf/1601.04937), p. 14, writes verbatim:

> "With no constraints, the convex hull ABCD has expected area √3 and expected perimeter
> (3+θ)√π = 6.4677562192310137839669010.... It is interesting that √3 is twice the expected
> area of a planar Gaussian triangle [17, 18, 19, 20], but the corresponding triangular
> perimeter 3√π is not so simply related to (3+θ)√π."

Two things follow. (a) The author treats the doubling as a numerical curiosity, not an
instance of a general identity — evidence that (I) is not in common circulation. (b) That
statement is about **Lebesgue area** for a Gaussian, whereas (I) is about the **`mu`-content**;
for a Gaussian these differ, so the two facts are logically distinct (both are true —
we checked both, `src/gaussian_identity_check.py`).

## 6. Search log (what was actually queried)

WebSearch: "Alikoski Sylvester four-point problem regular polygon exact formula expected area
random triangle"; "Marckert convex position five points disk probability p_5 explicit value";
"expected area convex hull four random points equals twice three random points convex body
identity Efron Renyi Sulanke"; "\"E[A_4] = 2 E[A_3]\" OR \"expected area of the convex hull of
four points is twice\""; "convex position five points regular polygon exact probability formula
n=5 Sylvester Buchta identity second moment triangle area"; "Buchta An identity relating moments
of functionals of convex hulls 2005"; "Buchta expected area convex hull n random points square
triangle exact 11/72 79/360 43/180"; "probability five points convex position formula second
moment triangle area covariance determinant"; "extremal conjecture convex position probability
triangle minimises ellipse maximises open problem n=5 proved"; "Around Sylvester's question in
the plane abstract new formula".
WebFetch (full text read): MathWorld *Polygon Triangle Picking*; arXiv:1402.3512 (Marckert,
PDF → text); arXiv:2411.08456 (Marckert–Morin, PDF → text); arXiv:1601.04937 (Finch, PDF → text);
arXiv:2401.16207 (abstract); labri.fr/perso/marckert/papers.html.
Blocked: Wiley (HTTP 402) and HAL (anti-bot) for the *published* Marckert–Rahmani 2021 text.
Resolved 2026-08-19 via the Semantic Scholar Graph API on the DOI, which exposed the otherwise
unadvertised preprint arXiv:1511.03658 (read in full). Lesson: when an author's page and the
citing papers give only a HAL/DOI link, query a metadata API on the DOI before giving up.

**Update 2026-08-19: that action is done.** Marckert & Rahmani was obtained via
arXiv:1511.03658 and read; identity (II) is not in it (see §4 and §5). The remaining
literature risk is narrow and stated: only the *published* 2021 text is unobtainable, and its
abstract describes the same "new formula" as the preprint.

**Next literature action, in priority order:** (1) library access to the published Mathematika
version, to close the last 5% of doubt on (II); (2) Buchta's variance papers
([*Exact formulae for variances of functionals of convex hulls*](https://www.plus.ac.at/wp-content/uploads/2021/02/Exact_formulae_for_variances_of_functionals_of_convex_hulls.pdf))
— the natural place for `E[A_3^3]` and `E[A_4^2]` of a general convex body, which is what
`n = 6` needs; (3) Alikoski 1939 itself (Ann. Acad. Sci. Fenn. 51 No. 7), still unobtained,
though its formula is now independently reproduced to 50 digits.


## 7. Addendum 2026-08-19 06:30 (local session)

* The n=5 theorem (Marckert-Rahmani 2021, section 4 above) now also has an independent short
  proof in this campaign via identity (II): `docs/N5_PROOF.md` (external model session; the
  operator relayed it; key lemma verified numerically here). The external session reports that
  Marckert-Rahmani's Theorems 1 and 3 give the bounds, equality cases, and the required
  Steiner/shaking monotonicities -- consistent with our reading of Marckert-Morin 2024 -- but
  we still have NOT read the paper ourselves, so the novelty status of identity (II) and of the
  section-lemma mechanism remains **UNRESOLVED**.
* Identities (I) and (II) also admit fully elementary proofs (deterministic four-point area
  identity; inclusion-exclusion on the number of non-vertices) recorded in `THEOREMS.md` section 8;
  these hold for every i.i.d. absolutely continuous planar law.


## 8. Marckert (-Rahmani) OBTAINED and read: arXiv:1511.03658 (2026-08-19 06:35, local session)

* The paper the overnight session could not get is on arXiv as **1511.03658** (v1, 11 Nov 2015,
  single author J.-F. Marckert; the Mathematika 67 (2021) 860-884 version adds S. Rahmani as
  coauthor). Found via the Semantic Scholar DOI record
  (https://api.semanticscholar.org/graph/v1/paper/DOI:10.1112/mtk.12110 -> externalIds.ArXiv);
  Wiley (402) and HAL (Anubis anti-bot) both still block. Text extraction saved as
  `docs/marckert_rahmani_1511.03658.txt` (22 pages, pymupdf).
* **What it does [verified, read in full]:** Theorem 1: `11/36 = Q^5_T <= Q^5_H <= Q^5_D = 1 - 305/(48 pi^2)`,
  equality only for triangles / ellipses. Theorem 3: for n = 4 and 5, Steiner symmetrisation
  (x-axis) does not decrease `Q^n_H` and shaking does not increase it. Method: condition on the
  abscissae (their law is invariant under both operations, eq. (10)); normalise by a
  verticality-preserving affine map so the extreme segments coincide (Def. 11); condition on the
  two extreme points and on which points lie above/below their line (Lemma 14); the two halves are
  "combs" and `<Cb[x_j, l_j]>` has an explicit rational recursive formula (Prop. 17, their "new
  formula for Q^n_H"); Prop. 13 then compares the resulting polynomial in the "symmetry defect"
  `beta` (`|q_j(beta)| <= p_j(lambda)`, `|beta_j| <= lambda_j`, cf. our `|d| <= r`), and for n = 5
  "remains only terms quadratic in beta, linear in lambda and in L" (p. 21) -- for n = 6 the
  polynomial has degree 2 and 4 in beta and they could not conclude.
* **What it does NOT contain [verified by full-text keyword scan and reading]:** no identity
  `P_5 = 1 - 10 E[A_3(1-A_3)]`, no `E[A_4] = 2 E[A_3]`, no covariance/inertia form, no
  Renyi-Sulanke/Efron/Buchta moment identity (Buchta is cited only for `Q^{n,m}`), no regular
  polygons, no half-disk. Its only "algebraic formula" in terms of areas is Blaschke's
  `Q^4 = 1 - 4 E(Area)`.
* **Novelty verdict, revised:** identity (II) and its covariance form are NOT in the arXiv version.
  The journal version (2021, with Rahmani, "some partial results" toward general n per its abstract)
  remains unread; presume it extends the comb formula rather than introducing moment identities,
  but this is UNVERIFIED. So (II), the regular-m-gon closed form (E4), and the fibrewise
  section-lemma proof in `docs/N5_PROOF.md` are "not found in the arXiv version + not found
  elsewhere" -- still not to be claimed as new without the journal text, but the presumption has
  moved from "probably known" to "possibly unstated".
* Relation between the two n = 5 proofs: same skeleton (Blaschke's method: fibres over abscissae,
  Steiner + shaking, Hausdorff limits). Marckert compares a rational/polynomial expression obtained
  from a 5-point comb recursion; the section-lemma proof uses the reduction to `E[A_3(1-A_3)]`
  (3-point fibres) and the flat-plateau density of `R = b_1U_1 - b_2U_2 + b_3U_3` to get the closed
  fibre formula `Phi(s) - Phi(t) = (s^2 - t^2)(V - b_2)/(4V^2 b_2)` -- the "terms quadratic in
  beta" of Marckert's p. 21, in one line, with the coefficient `(V - b_2)/(4V^2 b_2)` explicit.


## 9. Prior-art gate CLEARED (2026-08-19): triangle-area moments etc. -- full verdicts in docs/lit_triangle_area_moments.md

An external literature check (the Reed/Philip/Buchta gate) settled the moment-value novelty:
* **E[A_3^3](triangle)=31/9000 and (square)=137/72000 are KNOWN** (Reed 1974, Pacific J. Math.
  54:183-198; tabulated in Beck 2024, arXiv:2412.07952, Tables 2.2/2.4), as are E[A_3^4]
  (1/900, 1/2400). **Retire novelty for these; cite Reed/Beck.** (MathWorld's inline "137/9000"
  for the square is a typo; our 137/72000 = its own closed form at n=3.)
* **DISK higher triangle-area moments are NOT FOUND** -- MathWorld "Disk Triangle Picking" says
  the area distribution is "apparently not known exactly." So E[A_3^2]=3/(32 pi^2),
  E[A_3^3]=1001/(6400 pi^4), E[A_3^4](disk) are genuine novelty candidates.
* **E[A_4^2] / E[A_4^2 & convex] (Group 2): NOT FOUND** -- machinery exists (Efron 1965,
  Buchta 2005/2013) but not these values. Residual risk: Philip's `area12.pdf` (404, unread).
* **E[A_4] = 2 E[A_3] (any law): treat as KNOWN** (Efron 1965 / Buchta 1990 distribution-
  independent identities); cite, do not claim.
* **Identity (II) covariance form: NOT FOUND** (confirms the earlier verdict).
* **P_5(regular m-gon) and P_6(pentagon/hexagon/octagon): NOT FOUND** -- Morin 2024
  (arXiv:2401.16207) says exact regular-polygon formulas "are rare"; only asymptotics exist.
So the honest novelty candidates are: disk higher moments, E[A_4^2], identity (II), and the
regular-m-gon P_5/P_6 -- all "not found," pending Philip area12.pdf / Finch / Mathai p.391.
