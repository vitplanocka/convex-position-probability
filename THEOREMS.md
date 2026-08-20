# THEOREMS — statements, proofs and machine checks

*Written 2026-08-19 by the overnight server session. Everything here is either (a) a
classical result with a citation, or (b) proved in full below and checked by at least two
independent computations. Novelty verdicts are in `LITERATURE.md`; this file is about
correctness, not priority.*

## 0. Setting and notation

Let `mu` be an absolutely continuous probability measure on the plane. Draw
`X_1, X_2, ... ` i.i.d. from `mu` and put

    A_k := mu( conv(X_1, ..., X_k) )          (the mu-CONTENT of the hull, not its area)
    N_n := # vertices of conv(X_1, ..., X_n)
    P_n := P( X_1, ..., X_n are in convex position ) = P(N_n = n).

When `mu` is the uniform law on a convex body `K`, `A_k = area(conv)/|K|`, and `P_n = P_K(n)`.
All the quantities are invariant under invertible affine maps (they only involve ratios of
areas and incidence), so "square" means every parallelogram and "disk" every ellipse.

Define the **half-plane mass**

    c := mu( { x : x strictly left of the directed line X_1 -> X_2 } ).

`c` is a random variable in `[0,1]` determined by the (ordered) pair `(X_1, X_2)`.

Throughout, `m_k := E[c^k]` and `d := c - 1/2`, `mu_j := E[d^j]`.

---

## 1. Three classical ingredients

**Fact A (Renyi–Sulanke edge count).** For `n >= 3`,

    E[N_n] = C(n,2) * E[ c^{n-2} + (1-c)^{n-2} ].

*Proof.* The hull is a polygon, so its number of vertices equals its number of edges. The
segment `X_i X_j` is an edge of `conv(X_1..X_n)` exactly when the remaining `n-2` points all
lie strictly on one side of the line `X_i X_j`. Given the line, the remaining points are
i.i.d. `mu`, so that conditional probability is `c^{n-2} + (1-c)^{n-2}`. Sum over the `C(n,2)`
pairs and use exchangeability. (Ties have probability 0 because `mu` is absolutely
continuous.) ∎
[Renyi & Sulanke, Z. Wahrsch. 2 (1963) 75–84.]

**Fact B (Efron's identity).** For `n >= 2`,  `E[A_{n-1}] = 1 - E[N_n]/n`.

*Proof.* Condition on the `n` points and ask for `P(X_n not in conv(X_1..X_{n-1}))`. On one
hand this is `1 - E[A_{n-1}]`. On the other hand, `X_n` fails to be in the hull of the other
`n-1` exactly when `X_n` is a vertex of `conv(X_1..X_n)`; by exchangeability the probability
that a *given* one of the `n` points is a vertex is `E[N_n]/n`. ∎
[Efron, Biometrika 52 (1965) 331–343.]

**Fact C (Efron–Buchta inclusion identity).** For `3 <= m <= n`,

    P(N_n = m) = C(n,m) * E[ 1{X_1..X_m in convex position} * A_m^{n-m} ],

and consequently, writing `R(m,j) := E[1{conv position} A_m^j]`,

    R(k,j) = E[A_k^j] - sum_{m=3}^{k-1} C(k,m) R(m, j+k-m),
    P_n = R(n,0) = 1 - sum_{m=3}^{n-1} C(n,m) R(m, n-m)
        = sum_{j=0}^{n-3} (-1)^j C(n,j) E[ A_{n-j}^j ].

*Proof.* `N_n = m` means: some `m`-subset is in convex position **and** the other `n-m` points
lie inside its hull. The subset is unique, so summing over subsets is an exact
decomposition; conditionally on the chosen `m` points, each of the remaining `n-m` lands in
the hull independently with probability `A_m`. The recursion is the same statement read as
`E[A_k^j] = sum_{m=3}^{k} C(k,m) R(m, j+k-m)` (classify `conv(X_1..X_k)` by its vertex set)
and solved downwards. ∎
[Efron 1965; Buchta, *An identity relating moments of functionals of convex hulls*,
Discrete Comput. Geom. 33 (2005) 125–142; "Efron–Buchta identities".]

Special cases used below:

    P_4 = 1 - 4 E[A_3],                    (Sylvester)
    P(N_5 = 3) = 10 E[A_3^2],
    P(N_5 = 4) = 5 E[A_4] - 20 E[A_3^2],
    P_5 = 1 - 5 E[A_4] + 10 E[A_3^2].      (*)

`src/route_moments.py::coefficients(n)` prints these rational coefficients for any `n`
(derived programmatically from Fact C, not hard-coded).

---

## 2. The symmetry of the half-plane mass

**Lemma 1.** `c` and `1 - c` have the same distribution. Equivalently `d = c - 1/2` is
symmetric about 0, so **all odd moments of `d` vanish**. In particular `E[c] = 1/2`.

*Proof.* Swapping `X_1` and `X_2` reverses the direction of the line, exchanges the two open
half-planes, and hence replaces `c` by `1-c`. Since `(X_1, X_2)` and `(X_2, X_1)` have the same
joint law, `c =_d 1-c`. ∎

**Corollary 1.1.** The odd moments `m_1, m_3, m_5, ...` are determined by the even ones:
`m_1 = 1/2`, `m_3 = (3 m_2 - 1/2)/2`, `m_5 = (1 - 5 m_2 + 5 m_4)/2`, and in general
`2 m_{2k+1} = sum_{j=0}^{2k} (-1)^j C(2k+1, j) m_j`.

*Proof.* Expand `m_k = E[(1-c)^k]` binomially and solve for `m_k`; for odd `k` the coefficient
of `m_k` on the right is `-1`, giving `2m_k = sum_{j<k} (-1)^j C(k,j) m_j`. (For even `k` the
`m_k` terms cancel and one gets instead a *constraint* on the lower moments, the first of
which is exactly `m_1 = 1/2`.) ∎

---

## 3. The master formula and identity (I)

**Theorem 1 (master formula).** For every absolutely continuous `mu` on the plane and every
`n >= 3`,

    E[A_{n-1}] = 1 - (n-1) * E[c^{n-2}] = 1 - (n-1) * E[(1/2 + d)^{n-2}].

*Proof.* Facts A and B give `E[A_{n-1}] = 1 - C(n,2) E[c^{n-2}+(1-c)^{n-2}]/n
= 1 - (n-1)/2 * E[c^{n-2}+(1-c)^{n-2}]`, and by Lemma 1 the bracket is `2 E[c^{n-2}]`. ∎

Expanding `(1/2+d)^{n-2}` and killing the odd moments of `d` (Lemma 1):

    E[A_3] =   1/4  -   3    mu_2
    E[A_4] =   1/2  -   6    mu_2
    E[A_5] =  11/16 - (15/2) mu_2 -  5      mu_4
    E[A_6] =  13/16 - (15/2) mu_2 - 15      mu_4
    E[A_7] =  57/64 - (105/16) mu_2 - (105/4) mu_4 -  7 mu_6
    E[A_8] =  15/16 - (21/4)  mu_2 -  35      mu_4 - 28 mu_6
    E[A_9] = 247/256 - (63/16) mu_2 - (315/8) mu_4 - 63 mu_6 - 9 mu_8

**A new free parameter appears only at every OTHER step**: `E[A_3]` introduces `mu_2`,
`E[A_4]` introduces nothing new, `E[A_5]` introduces `mu_4`, `E[A_6]` introduces nothing new,
and so on. Hence:

**Theorem 2 (identity (I) and its family).** For every absolutely continuous `mu` on the plane

    E[A_4] = 2 E[A_3]                                            (I)
    E[A_6] = 3 E[A_5] - 5 E[A_3]                                 (I')
    E[A_8] = 4 E[A_7] - 14 E[A_5] + 28 E[A_3]                    (I'')

and in general every `E[A_{2j}]` is an explicit rational combination of
`E[A_3], E[A_5], ..., E[A_{2j-1}]`.

*Proof.* From the displayed expansions: `E[A_4] = 1/2 - 6 mu_2 = 2(1/4 - 3 mu_2) = 2 E[A_3]`;
`3 E[A_5] - 5 E[A_3] = 33/16 - (45/2) mu_2 - 15 mu_4 - 5/4 + 15 mu_2 = 13/16 - (15/2) mu_2
- 15 mu_4 = E[A_6]`; similarly for (I''). In general `E[A_{2j}] = 1 - 2j*E[c^{2j-1}]` and by
Corollary 1.1 `E[c^{2j-1}]` is a rational combination of `m_2, m_4, ..., m_{2j-2}`, which by
Theorem 1 are affine functions of `E[A_3], E[A_5], ..., E[A_{2j-1}]`. ∎

**Remark (what (I) does NOT say).** `A_k` is the `mu`-*content*, not the Lebesgue area. For a
Gaussian `mu` the two differ, and the Lebesgue-area doubling `E[|conv(4 Gaussian pts)|] =
sqrt(3) = 2 E[|conv(3 Gaussian pts)|]` observed by Finch (arXiv:1601.04937, p. 14) is a
*separate* statement — both hold, but Theorem 2 does not imply it. Both were checked
numerically here (`src/gaussian_identity_check.py`).

**Machine checks of Theorem 2** (`src/exact_rational_check.py`, exact `Fraction` arithmetic,
route P of §6 — no floating point anywhere):

| body | E[A_3] | E[A_4] | E[A_5] | E[A_6] | (I) | (I') | (I'') |
|---|---|---|---|---|---|---|---|
| triangle | 1/12 | 1/6 | 43/180 | 3/10 | 0 | 0 | 0 |
| square | 11/144 | 11/72 | 79/360 | 199/720 | 0 | 0 | 0 |
| aff.-regular hexagon | 289/3888 | 289/1944 | 149347/699840 | 62647/233280 | 0 | 0 | 0 |
| random rational pentagon | 68598907/915546564 | (=2x) | — | — | 0 | — | — |
| random rational heptagon | 17956351/241813452 | (=2x) | — | — | 0 | — | — |

symbolically on the true regular `m`-gons for `m = 3,4,5,6` (`src/mgon_exact.py exact`);
at 50 decimal digits for `m = 3..24, 30, 40, 60, 100` (residual `< 1e-49`); by Richardson-
extrapolated route P for the disk, a 3:1 ellipse, a half-disk and a stadium (ratio
`= 2.00000000000000`); and by Monte Carlo for the Gaussian, a skewed anisotropic Gaussian
mixture and a heavy-tailed `t_3` law (`|z| < 1.7`), confirming that convexity and uniformity
are genuinely not needed.

---

## 4. Identity (II): `P_5` from the first two moments of the triangle area

**Lemma 2.** For three i.i.d. points with covariance matrix `Sigma` in the plane,
`E[det(X_2-X_1, X_3-X_1)^2] = 6 det Sigma`.

*Proof.* Translation does not change `det`, so assume `E[X] = 0`. Writing `u_ij := X_i x X_j`
(scalar cross product), `det(X_2-X_1, X_3-X_1) = u_12 + u_23 + u_31`. For each term,
`E[u_12^2] = E[(x_1 y_2 - y_1 x_2)^2] = E[x_1^2]E[y_2^2] - 2E[x_1y_1]E[x_2y_2] + E[y_1^2]E[x_2^2]
= 2(Sigma_11 Sigma_22 - Sigma_12^2) = 2 det Sigma` by independence. Every cross term such as
`E[u_12 u_23]` contains exactly one factor that is linear in the independent, mean-zero
variable `X_3` (resp. `X_1`), hence vanishes. Three terms of `2 det Sigma` give `6 det Sigma`. ∎

**Corollary 2.1.** For the uniform law on a convex body `K`, since `A_3 = |det|/(2|K|)`,

    E[A_3^2] = (3/2) det(Sigma_K) / |K|^2 .

(Checks: square `1/96`, triangle `1/72`, disk `3/(32 pi^2)` — all reproduced exactly by
`polygon_exact.cov_det_area` and by Monte Carlo.)

**Theorem 3 (identity (II)).** For every absolutely continuous `mu` on the plane,

    P_5 = 1 - 10 E[A_3] + 10 E[A_3^2] = 1 - 10 E[ A_3 (1 - A_3) ],

and for the uniform law on a convex body `K`, equivalently

    P_K(5) = (5/2) P_K(4) - 3/2 + 15 det(Sigma_K)/|K|^2 .

*Proof.* Fact C gives `P_5 = 1 - 5 E[A_4] + 10 E[A_3^2]` (equation (*)). Substitute
`E[A_4] = 2 E[A_3]` from Theorem 2. For the second form use `P_4 = 1 - 4E[A_3]`, i.e.
`E[A_3] = (1-P_4)/4`, and Corollary 2.1. ∎

**Reading.** `P_5` depends on the body only through the **first two moments of the random
triangle's area**. The four-point problem controls the first moment; the second moment is a
pure covariance (an affine-invariant "isotropic constant"). No new integral is needed
beyond Sylvester's.

Machine checks: the exact `Fraction` table above returns `11/36` (triangle), `49/144`
(square, parallelogram), `1373/3888` (aff.-regular hexagon) — Valtr's values; the disk gives
`1 - 305/(48 pi^2)`; and the direct Monte-Carlo tester agrees on nine bodies
(`results/A3_identity_mc.json`).

---

## 5. Consequences

**Corollary 5.1 (disk).** `Sigma = (1/4) I`, `|K| = pi`, `E[A_3] = 35/(48 pi^2)`, so

    E[A_3^2] = 3/(32 pi^2),   E[A_4] = 35/(24 pi^2),   P_disk(5) = 1 - 305/(48 pi^2)
                                                                 = 0.3561883122726454...

This value is **not new**: it is in Marckert 2017 (arXiv:1402.3512), Table (7),
`1 - P^D_5 = 305/(48 pi^2)`. See `LITERATURE.md`.

**Corollary 5.2 (regular `m`-gon).** With `w = 2 pi / m`, Alikoski's 1939 four-point result
`E[A_3] = (9 cos^2 w + 52 cos w + 44)/(36 m^2 sin^2 w)` and the elementary
`E[A_3^2] = (2 + cos w)^2 / (24 m^2 sin^2 w)` (Corollary 2.1, using
`Sigma = ((2+cos w)/12) I` and `|K| = (m/2) sin w`) give

    P_5(regular m-gon) = 1 - 5 (15 cos^2 w + 92 cos w + 76) / (36 m^2 sin^2 w).

Checks built into `src/mgon_table.py`: `m = 3` gives `11/36`, `m = 4` gives `49/144` (both
Valtr), and `m -> infinity` gives `1 - 305/(48 pi^2)` (the disk). See `docs/N5_LANDSCAPE.md`.

**Corollary 5.3 (`E[N_6] = 30 E[c^4]`, and what fails at `n = 6`).** By Lemma 1,
`T_4 = E[c^4 + (1-c)^4] = 2 E[c^4]`, so Fact A gives

    E[N_6] = C(6,2) T_4 = 30 E[c^4],   and   E[A_5] = 1 - E[N_6]/6 = 1 - 5 E[c^4].

(For the disk this is a one-dimensional integral; `src/disk_EA5_exact.py` evaluates it
symbolically and gets `E[N_6](disk) = 6 - 175/(12 pi^2) + 23023/(1152 pi^4)`, hence
`E[A_5](disk) = 7(2400 pi^2 - 3289)/(6912 pi^4) = 0.212072074036700971...`, which agrees
with the completely independent route-P Richardson value to `2.9e-17`, and agrees *exactly*
with `E[N_6]` recomputed from Marckert 2017's `P^D_{6,m}` table.)

Now `P_6 = 1 - 6 E[A_5] + 15 E[A_4^2] - 20 E[A_3^3]`. The first term needs `mu_4`, a
genuinely new parameter not determined by `E[A_3]`; and `E[A_4^2]`, `E[A_3^3]` are
second/third moments of hull areas, which the edge-count route does not reach at all. So the
`n = 5` collapse is special: it happens because `E[A_4]` is the one even-index hull moment
that introduces no new parameter. (For the disk all three are nevertheless known — see
`PROGRESS.md` 2026-08-19 00:19: `E[A_3^3] = 1001/(6400 pi^4)` and
`E[A_4^2] = (2400 pi^2 + 31031)/(19200 pi^4)` follow from Marckert's `P^D_{6,3}` and
`P^D_{6,4}` through Fact C, and the three assemble to his `P_6` exactly.)

---

## 5b. A three-dimensional analogue

The planar proof of (I) used `#vertices = #edges`, which is false in `R^3`. But for a
**simplicial** 3-polytope Euler's relation `V - E + F = 2` together with `3F = 2E` (every
facet a triangle, every edge in two facets) gives exactly

    V = 2 + F/2,

and points in general position give a simplicial hull. The Renyi–Sulanke facet count in
`R^3` reads `E[F_n] = C(n,3) E[c^{n-3} + (1-c)^{n-3}]`, where now `c` is the `mu`-mass on one
side of the plane through *three* of the points; `c` and `1-c` are again identically
distributed (a transposition of the three points reverses the orientation, and the law is
exchangeable), so Lemma 1 and Corollary 1.1 carry over verbatim with `m_k = E[c^k]`. Hence
`E[F_n] = 2 C(n,3) m_{n-3}`, `E[N_n] = 2 + C(n,3) m_{n-3}`, and Efron gives

    E[A_{n-1}] = 1 - 2/n - C(n,3) m_{n-3} / n :
      E[A_3] = 1/2 - m_1   = 0        (correct — three points span no volume)
      E[A_4] = 3/5 - 2 m_2
      E[A_5] = 2/3 - (10/3) m_3 = 3/2 - 5 m_2
      E[A_6] = 5/7 - 5 m_4
      E[A_7] = 3/4 - 7 m_5 = -11/4 + (35/2) m_2 - (35/2) m_4

**Theorem 4 (identity (I_3)).** For every absolutely continuous probability measure on `R^3`,

    E[A_5] = (5/2) E[A_4] ,   and   E[A_7] = (7/2) E[A_6] - (35/4) E[A_4].

*Proof.* Substitute `m_3 = (3 m_2 - 1/2)/2` (Corollary 1.1) into the `E[A_5]` line above to
get `3/2 - 5 m_2 = (5/2)(3/5 - 2 m_2) = (5/2) E[A_4]`; similarly `m_5 = (1 - 5m_2 + 5m_4)/2`
gives the second relation. ∎

The same "one new parameter every other step" pattern therefore holds in `R^3`, with the
constants shifted by the `V = 2 + F/2` bookkeeping. Combined with Efron
(`P_5^{(3)} = 1 - 5 E[A_4]`), (I_3) says `E[A_5]` in three dimensions is determined by the
3-D five-point Sylvester problem alone. Anchors and consequences: `E[A_4](ball) = 9/715` gives
**`E[A_5](ball) = 9/286 = 0.0314685314...`** and `E[A_4](cube) = 3977/216000 - pi^2/2160`
(Zinani 2003) gives **`E[A_5](cube) = 3977/86400 - pi^2/864 = 0.0346069393...`** — i.e. the
expected volume of the hull of *five* uniform points comes free from the *four*-point value. Numerical test: `src/dim3_identity.py`
(`results/dim3_identity.json`), which computes `N_5` and `N_6` by exact orientation
predicates — a point is a non-vertex iff it lies in the hull of the others, which in `R^3`
means inside `conv` of four of them (Caratheodory) — so no hull construction and no
floating-point volume enters.

---

## 6. Route P: exact Sylvester quantities for an arbitrary convex polygon

**Proposition 6.1.** Let `K` be a convex polygon with vertices `V_0, ..., V_{m-1}` (CCW),
edges `E_i = [V_i, V_{i+1}]`, `d_i = V_{i+1} - V_i`. For each unordered pair `i < j`
parametrise a line crossing `E_i` and `E_j` by its two boundary points
`P(u) = V_i + u d_i`, `Q(v) = V_j + v d_j`, `(u,v) in [0,1]^2`. Then for every `k >= 0`

    J_k := \int_{lines} L^3 [ c^k + (1-c)^k ] dG
         = sum_{i<j} \int_0^1 \int_0^1 W(u,v) [ c(u,v)^k + (1-c(u,v))^k ] du dv,
    W(u,v) = - ( d_i x (Q - V_i) ) * ( d_j x (V_j - P) )  >= 0,
    c(u,v) = shoelace( [P, V_{i+1}, ..., V_j, Q] ) / |K|,

where `L` is the chord length, `dG = dp dtheta` the rigid-motion-invariant line measure and
`c` the fraction of `|K|` cut off. The integrand is a **polynomial** of degree `<= 2k+2`
in `(u, v)`, so `J_k` is elementary. Moreover

    T_k := J_k / (3 |K|^2) = E[ c^k + (1-c)^k ],   E[N_n] = C(n,2) T_{n-2}.

*Proof.* Two i.i.d. uniform points determine a line; the planar Blaschke–Petkantschin /
Cauchy formula gives `dx_1 dx_2 = |t_1 - t_2| dt_1 dt_2 dG`, and integrating the two
positions along a chord of length `L` gives `\int\int |t_1-t_2| dt_1 dt_2 = L^3/3`. Hence for
any `f`, `E[f(c)] = (1/(3|K|^2)) \int L^3 f(c) dG`; taking `f = 1` and using
`\int L^3 dG = 3|K|^2` fixes the normalisation. Santalo's boundary parametrisation of the
line measure reads `dp dtheta = sin(phi_i) sin(phi_j) ds_i ds_j / L` with `phi` the angles
between the chord and the boundary; with `w = (Q-P)/L`, `sin(phi_i) ds_i = |e_i x w| |d_i| du
= |d_i x (Q-P)| du / L`, and likewise for `j`, so `dG = |d_i x (Q-P)||d_j x (Q-P)| du dv / L^3`
and the `L^3` cancels against the `L^3` in `J_k`. Finally `d_i x (Q-P) = d_i x (Q - V_i)`
depends only on `v` (because `d_i x (P - V_i) = 0`) and has a constant sign over the cell by
convexity; likewise `d_j x (Q-P) = d_j x (V_j - P)` depends only on `u`. Every line meeting
the interior crosses exactly two edges, so the cells tile the line space. ∎

**Built-in exactness checks** (both must hold identically, and do):
`T_0 = 2` (equivalently `\int L^3 dG = 3|K|^2`) and `E[N_3] = 3 T_1 = 3`
(three points always span a triangle). Independent validation: route P reproduces
Alikoski's regular-`m`-gon values `1/12, 11/144, (9+2 sqrt 5)/180, 289/3888,
(97+52 sqrt 2)/2304` for `m = 3,4,5,6,8` and agrees with Alikoski to 50 digits for
`m = 3..24, 30, 40, 60, 100`; on polygonal approximations plus Richardson extrapolation
(error `~ C/m^4`, verified) it returns the disk's `35/(48 pi^2)` to `6e-17`.

---

## 7. Summary of what is proved vs. assumed

| statement | status |
|---|---|
| Facts A, B, C | classical, cited, and each re-derived above |
| Lemma 1 (`c =_d 1-c`) | proved (one line) |
| Theorem 1 (master formula) | proved from A + B + Lemma 1 |
| Theorem 2 = identity (I) and family | proved; checked exactly on 6 polygons, symbolically on 4 regular `m`-gons, to 50 digits on 26 more, and by MC on 3 non-uniform laws |
| Lemma 2 / Corollary 2.1 | proved; checked exactly |
| Theorem 3 = identity (II) | proved from Fact C + Theorem 2 |
| Corollary 5.1 (`P_5` disk) | proved; **known** (Marckert 2017) |
| Corollary 5.2 (`P_5` regular `m`-gon) | proved; novelty: see `LITERATURE.md` |
| Proposition 6.1 (route P) | proved; two identical-by-construction checks pass |
| Theorem 4 (n=5 extremal) | **known** (Marckert-Rahmani 2021); short independent proof via (II) in `docs/N5_PROOF.md`, section lemma verified numerically to 9 digits |

---

## 8. Theorem 4 (n = 5 extremal; Marckert-Rahmani 2021) -- a short proof via identity (II)

**Theorem 4.** For every planar convex body `K`, `11/36 <= P_5(K) <= 1 - 305/(48 pi^2)`, with
equality on the left only for triangles and on the right only for ellipses.

Full proof, verification record and attribution: [`docs/N5_PROOF.md`](docs/N5_PROOF.md)
(external model session, relayed by the operator 2026-08-19; re-derived and numerically
verified here). Mechanism: with `F = E[A_3(1-A_3)]` and vertical sections, on each ordered
triple of abscissae the noise `R = b_1 U_1 - b_2 U_2 + b_3 U_3` has a *flat* density `1/(2 b_2)`
on `[-r, r]` (r = b_2 - b_1 - b_3 >= 0 by convexity), so the fibre contribution satisfies
`Phi(s) - Phi(t) = (s^2 - t^2)(V - b_2)/(4 V^2 b_2)` with `b_2 < V`; Steiner symmetrisation
sends the centre second-difference `d -> 0` and shaking sends `d -> -r`, giving exact
nonnegative gap formulas (3.1)/(3.2) there. Fibre by fibre `Delta E[A^2]/Delta E[A] = b_2/V < 1`,
which is why `F` inherits Blaschke's monotonicity even though `E[A_3]` and `E[A_3^2]` move the
same way (this answers the doubt recorded in `docs/N5_LANDSCAPE.md` section 6).

The same source gives **elementary re-proofs of Theorems 2 and 3** that need none of Facts A-C:
`2|conv(z_1..z_4)| = sum_i |conv{z_j: j != i}|` for any four points in general position (so
`E[A_4] = 2E[A_3]` for every i.i.d. law, deterministically), and
`1_{S=0} = 1 - S + C(S,2)` for `S` = number of non-vertices among five points (`S in {0,1,2}`),
which yields `P_5 = 1 - 5E[A_4] + 10E[A_3^2]` directly.
