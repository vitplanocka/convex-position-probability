# Prior-art / novelty check: triangle-area moments, hull-area moments, convex-position probabilities

*External literature check delivered to the operator 2026-08-19 and integrated by the local
session. This is the "Reed/Philip/Buchta gate" the campaign flagged. Bottom line below; full
per-group findings follow. One residual source (Philip's `area12.pdf`) 404s and could not be read.*

## Bottom line (novelty verdicts)

| result | verdict | source |
|---|---|---|
| E[A_3^3](triangle) = 31/9000 | **KNOWN** | Reed 1974; Beck 2024 (arXiv:2412.07952) Table 2.2 |
| E[A_3^3](square) = 137/72000 | **KNOWN** | Reed 1974; Beck 2024 Table 2.4; = MathWorld's closed form at n=3 |
| E[A_3^4](triangle)=1/900, (square)=1/2400 | **KNOWN** | Beck 2024 Tables 2.2/2.4 |
| E[A_3^2](disk)=3/(32 pi^2), E[A_3^3](disk)=1001/(6400 pi^4), E[A_3^4](disk) | **NOT FOUND** | MathWorld "Disk Triangle Picking": *"The distribution of areas ... is apparently not known exactly."* |
| E[A_4^2] and E[A_4^2 & convex] (triangle/square/disk) | **NOT FOUND** | machinery exists (Efron 1965, Buchta 2005/2013), values not located |
| E[A_4] = 2 E[A_3] as a general (any-law) identity | **KNOWN-ish** | consequence of Efron 1965 / Buchta 1990 distribution-independent identities; cite, do not claim |
| P_5 = 1 - 10 E[A_3(1-A_3)] = (5/2)P_4 - 3/2 + 15 detSigma/V^2 | **NOT FOUND** | not in Marckert-Rahmani (checked) nor elsewhere; covariance term unlocated |
| P_5(regular m-gon) closed form; P_6(pentagon/hexagon/octagon) | **NOT FOUND** | Valtr = triangle/parallelogram only; Marckert = disk; Morin 2024: *"exact formulas are rare"* (asymptotics only) |

**So: retire novelty for the triangle/square E[A_3^3] and E[A_3^4] (cite Reed 1974 / Beck 2024).
The genuine novelty candidates are (i) the DISK higher triangle-area moments, (ii) the 4-point
hull-area second moment E[A_4^2] / E[A_4^2&convex], (iii) the covariance-form P_5 identity, and
(iv) the exact regular-m-gon P_5 and P_6.** "Not found" = "not located in accessible literature,"
not "provably absent."

## Per-group detail

### Group 1 -- random-triangle area moments
- Triangle E[A_3^3] = 31/9000 and E[A_3^4] = 1/900: **KNOWN**. Beck (arXiv:2412.07952, "On Random
  Simplex Picking Beyond the Blaschke Problem," 2024) Table 2.2 lists v_2^(k)(T_2) = 1, 1/12, 1/72,
  31/9000, 1/900, ...; Beck credits Reed 1974 / Mathai p.391 / Maesumi. Reed, "Random points in a
  simplex," Pacific J. Math. 54 (1974) 183-198: abstract gives "explicit expressions ... for the
  moments of the area of the triangle formed by three points chosen at random in a triangle or
  parallelogram."
- Square E[A_3^3] = 137/72000 and E[A_3^4] = 1/2400: **KNOWN**. Beck Table 2.4. Also equals
  MathWorld's closed form mu'_n = 3*2^(3-n)[(n+2)H_{n+1}+1]/[(n+1)(n+2)^3(n+3)^2] at n=3
  (verified: H_4=25/12 -> 137/72000). NOTE: MathWorld's *inline* list misprints this as 137/9000;
  our 137/72000 is correct and matches Beck + the formula.
- Disk E[A_3^2]=3/(32 pi^2), E[A_3^3]=1001/(6400 pi^4), E[A_3^4]: **NOT FOUND**. MathWorld "Disk
  Triangle Picking" gives only the mean 35/(48 pi^2) and states the area distribution is "apparently
  not known exactly." Searched MathWorld, Philip's KTH reports, Finch, OEIS, arXiv -- absent.

### Group 2 -- 4-point hull-area second moment: NOT FOUND
The framework is Efron 1965 (Biometrika 52:331-343, area<->vertex-count identity) extended to
higher moments by Buchta 2005 (DCG 33:125-142). Buchta 2013 ("Exact formulae for variances of
functionals of convex hulls," Adv. Appl. Prob. 45:917-924) gives exact variances but for the
asymptotic corner functionals, not E[A_4^2] for n=4 in triangle/square/disk. Philip's "The Area of
a Random Convex Polygon" (TRITA-MAT-04-MA-07; `area12.pdf`) is the most likely place a square
E[A_4^2] could pre-exist -- **it 404s and could not be read (residual risk).** Otherwise variance of
the n-point hull area is treated only asymptotically (Pardon 2011, Ann. Probab. 39:881-903).

### Group 3 -- identities
- E[A_4] = 2 E[A_3] (any law): **KNOWN-ish** -- a distribution-independent linear identity of the
  Efron 1965 / Buchta 1990 ("Distribution-independent properties of the convex hull of random
  points," J. Theoret. Probab. 3:387-393) type (valid for n-d >= 2 and even; d=2,n=4 qualifies).
  The exact coefficient-2 statement was not extracted verbatim, but treat as known and cite.
- P_5 covariance identity: **NOT FOUND** (checked Marckert-Rahmani, Sylvester-problem literature,
  Barany 1999, Valtr, Marckert 2017, arXiv). The 15 detSigma/V^2 term unlocated anywhere.

### Group 4 -- regular-polygon P_5, P_6: NOT FOUND (exact)
Valtr 1995/96 = parallelogram (49/144) and triangle (11/36) only; Marckert 2017 = disk (recursive);
Morin 2024 (arXiv:2401.16207, 2410.11706) = regular-kappa-gon ASYMPTOTICS only, explicitly "exact
formulas are rare." No exact P_5 or P_6 for the regular pentagon/hexagon/octagon/general m-gon
found. Consistency check passed: our hexagon E[A_3] = 289/3888 matches MathWorld's Alikoski-based
"Hexagon Triangle Picking" value, and hexagon P_4 = 683/972.

## New references worth having
- M. Beck, "On Random Simplex Picking Beyond the Blaschke Problem," arXiv:2412.07952 (2024) -- has
  the simplex area-moment tables (Reed's moments, and beyond); directly relevant.
- L. Morin, "Probability that n points are in convex position in a regular kappa-gon: asymptotic
  results," arXiv:2401.16207 (2024); "... general convex polygon," arXiv:2410.11706 (2024).
- W. J. Reed, Pacific J. Math. 54 (1974) 183-198. B. Efron, Biometrika 52 (1965) 331-343.
- C. Buchta: DCG 33 (2005) 125-142; J. Theoret. Probab. 3 (1990) 387-393; Adv. Appl. Prob. 45
  (2013) 917-924. J. Pardon, Ann. Probab. 39 (2011) 881-903.
- Philip KTH reports (johanph): `area12.pdf` (404), `squaref.pdf`, `hexagon.pdf`, `area21.pdf`.

## Residual risk (before any "new" claim in a write-up)
1. Philip `area12.pdf` / `squaref.pdf` (square hull-area / triangle-area distributions) -- the most
   likely prior source for the square E[A_4^2]; `area12.pdf` 404s.
2. Finch, "Mathematical Constants II" (2019) random-triangle essays.
3. Mathai, "An Introduction to Geometrical Probability" (1999) p. 391 (general simplex moments).
A MathSciNet/zbMATH full-text search would further harden the "not found" verdicts.
