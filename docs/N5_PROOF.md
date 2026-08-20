# A short proof of the n = 5 extremal theorem via identity (II)

*Integrated 2026-08-19 06:30. The argument below was produced by an external model session
(prompted with `docs/PROBLEM_B.md`) and relayed by the operator; it was then re-derived line by
line and its key lemma verified numerically here (`src/section_lemma_check.py`, results at the
end). The THEOREM itself is not new: it is Marckert & Rahmani, "Around Sylvester's question in
the plane", Mathematika 67 (2021) 860-884 (Theorems 1 and 3 there, by a different route -- a
general formula for Q_H^n and a polynomial comparison under symmetrisation). What is recorded
here is a considerably shorter proof made possible by the reduction P_5 = 1 - 10 E[A_3(1-A_3)],
together with two elementary re-proofs of identities (I) and (II) that need none of
Renyi-Sulanke / Efron / Buchta.*

Notation: K a planar convex body, V = |K|; X_1, X_2, ... i.i.d. uniform in K;
A = A_3 = area(conv(X_1,X_2,X_3))/V;  P_n(K) = P(X_1..X_n in convex position);
F(K) = E[A(1-A)] = E[A] - E[A^2].

---

## Theorem A (n = 5 extremal; Marckert-Rahmani 2021)

For every planar convex body K,

    11/36  <=  P_5(K)  <=  1 - 305/(48 pi^2),

with equality on the left iff K is a triangle and on the right iff K is an ellipse.
Equivalently  61/(96 pi^2) <= F(K) <= 5/72  with the same equality cases.

The proof has four steps: (1) the identity P_5 = 1 - 10 F; (2) a fibrewise "section lemma"
giving EXACT gap formulas for Steiner symmetrisation and shaking; (3) monotonicity of F under
both operations; (4) Blaschke's convergence argument plus his uniqueness theorem.

---

## 1. Direct derivation of P_5 = 1 - 10 E[A(1-A)]  (no edge counts needed)

Sample five points; almost surely they are in general position. For i = 1..5 let
C_i = {X_i in conv{X_j : j != i}}  ("X_i is not a hull vertex") and S = sum_i 1_{C_i}.
The hull has 5, 4 or 3 vertices, so S in {0,1,2} a.s., and on these values

    1_{S=0} = 1 - S + C(S,2)         (check: S=0 -> 1; S=1 -> 0; S=2 -> 1-2+1 = 0).

Taking expectations,

    P_5 = 1 - sum_i P(C_i) + sum_{i<j} P(C_i and C_j).                               (1.1)

**One-point terms.** Conditioning on the other four points, P(C_i) = E[H_4] where
H_4 = area(conv of four points)/V. For ANY four points z_1..z_4 in general position,

    2 |conv(z_1..z_4)| = sum_{i=1}^4 |conv{z_j : j != i}|                            (1.2)

(convex quadrilateral: the four triangles pair into the two diagonal decompositions;
one point inside the triangle of the other three: the three small triangles partition the
outer one, plus the outer one itself). Taking expectations, 2 E[H_4] = 4 E[A], i.e.

    E[A_4] = 2 E[A_3]      -- identity (I), now for EVERY i.i.d. law, deterministically.  (1.3)

**Two-point terms.** C_i and C_j together force S = 2, i.e. the hull is the triangle T of
the other three points and both X_i, X_j lie in T; conversely that implies both. Conditioning
on the three others, P(C_i and C_j) = E[A^2].

Substituting into (1.1):  P_5 = 1 - 5 E[H_4] + 10 E[A^2] = 1 - 10 E[A] + 10 E[A^2], i.e.

    P_5 = 1 - 10 E[A(1-A)].                                                          (1.4)

**Covariance form.** With Y_i = X_i - E[X] and D = det(X_2-X_1, X_3-X_1) = det(Y_2,Y_3) +
det(Y_3,Y_1) + det(Y_1,Y_2): for independent centred Y, Z with covariance Sigma,
E[det(Y,Z)^2] = 2 det Sigma, and the cross terms vanish by independence and centring, so
E[D^2] = 6 det Sigma. Since A = |D|/(2V),  E[A^2] = (3/2) det Sigma / V^2, and with
P_4 = 1 - 4E[A]:

    P_5(K) = (5/2) P_4(K) - 3/2 + 15 det(Sigma_K)/|K|^2 .                             (1.5)

---

## 2. The section lemma

Fix a direction and write K in vertical-section form
K = {(x,y): a <= x <= b, l(x) <= y <= u(x)}, u concave, l convex. Put
c(x) = (u+l)/2 (centre), h(x) = (u-l)/2 (half-length), so the section at x is
[c-h, c+h] and V = 2 int_a^b h. A uniform point has x-density 2h(x)/V and, given x,
y = c(x) + h(x) U with U ~ Unif[-1,1].

Order the three abscissae x_1 < x_2 < x_3, put p = x_2 - x_1, q = x_3 - x_2,
c_i = c(x_i), h_i = h(x_i), and define

    d  = q c_1 - (p+q) c_2 + p c_3           (second difference of the centre line)
    b_1 = q h_1,  b_2 = (p+q) h_2,  b_3 = p h_3,     r = b_2 - b_1 - b_3.

**2.1 Convexity restrictions.**  x_2 = (q x_1 + p x_3)/(p+q) is a convex combination, so
concavity of u = c+h gives  q u_1 - (p+q) u_2 + p u_3 <= 0, i.e.  d - r <= 0;
convexity of l = c-h gives  q l_1 - (p+q) l_2 + p l_3 >= 0, i.e.  d + r >= 0.  Hence

    |d| <= r,   in particular  b_2 >= b_1 + b_3.                                    (2.1)

Moreover h is concave and nonnegative, so int_a^{x_2} h >= (x_2-a) h_2/2 and
int_{x_2}^b h >= (b-x_2) h_2/2, whence V >= (b-a) h_2 > (x_3-x_1) h_2 = b_2 a.s.:

    b_2 < V.                                                                         (2.2)

**2.2 The signed area.**  The doubled signed area of the ordered triple is
D = q y_1 - (p+q) y_2 + p y_3 = d + R  with  R = b_1 U_1 - b_2 U_2 + b_3 U_3, U_i i.i.d.
Unif[-1,1]. R is symmetric. Put W = b_1 U_1 + b_3 U_3, |W| <= b_1 + b_3; given W = w,
R = w - b_2 U_2 is uniform on [w - b_2, w + b_2], which contains [-r, r] for every w. Hence

    f_R(t) = 1/(2 b_2)   for all |t| <= r     (an exactly flat plateau),               (2.3)

and every admissible value -d lies on the plateau because |d| <= r. Let m(t) = E|R + t|.
For |t| <= r, m'(t) = 2 F_R(t) - 1 = t/b_2, so

    m(t) - m(0) = t^2 / (2 b_2),        and        E[(R+t)^2] - E[R^2] = t^2.          (2.4)

Define the fibre contribution to F,  Phi(t) = E[ |R+t|/(2V) - (R+t)^2/(4V^2) ].  Then for
s, t in [-r, r]:

    Phi(s) - Phi(t) = (s^2 - t^2) (V - b_2) / (4 V^2 b_2)   >= 0 whenever |s| >= |t|,   (2.5)

because b_2 < V. Equivalently, with M(t) = E|R+t|/(2V) the fibre contribution to E[A]:

    Phi(s) - Phi(t) = (1 - b_2/V) (M(s) - M(t)).                                     (2.6)

The joint density of the ordered abscissae is 6 prod (2h_i/V) = 48 h_1 h_2 h_3 / V^3.

---

## 3. Monotonicity

**Steiner symmetrisation** S K replaces each section by [-h(x), h(x)]: h unchanged, c -> 0,
so d -> 0 and nothing else in the fibre changes. By (2.5) with (s,t) = (d,0):

    F(K) - F(SK) = (12/V^5) int_{a<x_1<x_2<x_3<b} h_1 h_2 h_3 (d^2 (V-b_2)/b_2) dx  >= 0.   (3.1)

**Shaking** Sh K replaces each section by [0, 2h(x)]: h unchanged, c -> h, so d -> -r. By
(2.5) with (s,t) = (-r, d) and |d| <= r:

    F(ShK) - F(K) = (12/V^5) int h_1 h_2 h_3 ((r^2 - d^2)(V-b_2)/b_2) dx  >= 0,          (3.2)

where r^2 - d^2 = (r-d)(r+d) is the product of the concavity deficit of u and the convexity
deficit of l on the triple. In terms of P_5 = 1 - 10F:  P_5(SK) >= P_5(K),  P_5(ShK) <= P_5(K).

---

## 4. Conclusion and equality cases

Classically there is a sequence of Steiner symmetrisations of K converging in Hausdorff
distance to a disk of area V, and a sequence of shakings converging to a triangle of area V
(Blaschke; for shaking see also Campi-Colesanti-Gronchi 1999). F is continuous under
Hausdorff convergence of convex bodies with nonempty interior (bounded integrands over K^3;
boundaries are null). Hence F(K) >= F(disk) and F(K) <= F(triangle), i.e. Theorem A.

Equality: by (2.6) each fibre's F-deficit is (1 - b_2/V) times its E[A]-deficit with a
strictly positive multiplier a.e., so if F(K) = F(disk) every Steiner deficit of F vanishes,
hence every deficit of E[A] vanishes, hence E[A](K) = E[A](disk) and K is an ellipse by the
equality case of Blaschke's theorem; likewise F(K) = F(triangle) forces K to be a triangle.

**The finer inequality behind the result.** Fibre by fibre, under either operation,

    Delta E[A^2] / Delta E[A] = b_2 / V < 1,

so the first triangle-area moment always changes faster than the second, and F = E[A] - E[A^2]
inherits Blaschke's monotonicity. More generally F_lambda = E[A] - lambda E[A^2] is
nonincreasing under symmetrisation and nondecreasing under shaking for every lambda <= 1;
lambda = 1 is exactly the coefficient the geometric bound b_2 <= V justifies.

Endpoints: F(triangle) = 1/12 - 1/72 = 5/72 (P_5 = 11/36); F(disk) = 35/(48 pi^2) - 3/(32 pi^2)
= 61/(96 pi^2) (P_5 = 1 - 305/(48 pi^2)); F(square) = 11/144 - 1/96 = 19/288 (P_5 = 49/144).

---

## 5. Verification record (this campaign, 2026-08-19)

`src/section_lemma_check.py`:

* (2.3)-(2.4) on one fibre (b = (0.3, 1, 0.25), r = 0.45; 4e6 samples): density 0.495-0.504
  vs claimed 0.500 across the plateau; m(t)-m(0) = 0.08005, 0.02002, 0, 0.01997, 0.08002 vs
  claimed 0.08, 0.02, 0, 0.02, 0.08.
* (3.1) on the half-disk {x^2+y^2 <= 1, y >= 0} with vertical sections (SK = the ellipse
  with semi-axes 1, 1/2, F(SK) = F(disk)): the deficit integral, Sobol QMC 2^22 points x 4
  seeds, gives 0.001511838 +- 5.5e-11.
* (3.2) on the disk (d = 0; ShK = a half-ellipse, affinely a half-disk): 0.001511838 +- 5.5e-11.
* Target: F(half-disk) - F(disk) with E[A_3](half-disk) = 0.076512497523552 (route P +
  Richardson, an unrelated method) and E[A_3^2](half-disk) = (3/2) det Sigma/V^2 exact
  = 0.010619491187: **0.001511837564**. Three-way agreement to 9 digits.
  Also P_5(half-disk) = 1 - 10 F = 0.341069936632, matching the overnight value
  0.341069936631669 (VERDICT.md, A3 table).
* Sanity in the disk integral: min r = 7.8e-13 (>= 0), max b_2/V = 0.636 (< 1).

Status: Theorem A is a published theorem (Marckert-Rahmani 2021). The proof above is an
independent, shorter route; identities (1.3)-(1.5) are proved here in full generality (any
i.i.d. absolutely continuous planar law) by elementary means. Whether (1.5) or the
section-lemma mechanism appears in Marckert-Rahmani remains UNRESOLVED until that paper is
read (see LITERATURE.md section 4).


---

## 6. Update 2026-08-19 06:40 -- Marckert's paper obtained (arXiv:1511.03658)

The arXiv version was read in full (see LITERATURE.md section 8). It uses the same Blaschke
skeleton (fibres over abscissae, Steiner symmetrisation and shaking, Hausdorff limits) but a
different engine: a normalised 5-point "comb" recursion (its Props. 16-17) whose polynomial in the
symmetry defect beta is compared for n = 4, 5 ("for n = 5 remains only terms quadratic in beta",
p. 21). It contains no moment identity, no `E[A_4] = 2E[A_3]`, and no covariance form. So the proof
above is genuinely a different (and shorter) route to the same theorem; the identity-based
reduction and the closed fibre formula are not in that text. The 2021 journal version with Rahmani
is still unread; novelty remains formally UNRESOLVED but is now "not found in the arXiv version".
