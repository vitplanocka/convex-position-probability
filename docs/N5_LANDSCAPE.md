# The n = 5 landscape

*Server session, 2026-08-19. Everything here is computed with route P
(`src/polygon_exact.py`), which gives `E[A_3]` — hence `P_4` and, via identity (II),
`P_5` — **exactly** for any convex polygon. No Monte Carlo enters this document except
as an independent cross-check.*

## 0. The one thing to read first

**The n = 5 extremal conjecture is not open.** Marckert & Rahmani, *Around Sylvester's
question in the plane*, Mathematika 67 (2021) 860–884, proved
`P_triangle(5) <= P_K(5) <= P_disk(5)` for every planar convex body `K`; `n >= 6` remains
open. (Quoted verbatim in Marckert & Morin, [arXiv:2411.08456](https://arxiv.org/pdf/2411.08456),
p. 4.) TASK.md section B3 asked us to test that conjecture and hunt for a proof; what
follows is therefore an **independent numerical verification of a published theorem**, plus
an exact table that we could not find anywhere in the literature.

## 1. Why n = 5 is easy to survey

By identity (II) (`THEOREMS.md`, Theorem 3),

    P_5(K) = 1 - 10 F(K),      F(K) := E[A_3] - E[A_3^2],

so the whole n = 5 problem is the two-sided inequality

    61/(96 pi^2) = 0.064381168773  <=  F(K)  <=  5/72 = 0.069444444444,

with equality for ellipses (left) and triangles (right). `E[A_3]` is Sylvester's four-point
functional; `E[A_3^2] = (3/2) det Sigma_K / |K|^2` is a pure covariance. Both are maximised
by triangles and minimised by ellipses (Blaschke 1917; the planar isotropic-constant
extremes), so `F` is a **difference of two same-direction functionals** and the inequality is
not a formal consequence of either — which is presumably why n = 5 needed a real proof.

## 2. Exact table for regular m-gons

`E[A_3]` from Alikoski 1939 — independently reproduced by route P to 50 digits, see
`LITERATURE.md` §2 — and `E[A_3^2] = (2 + cos w)^2/(24 m^2 sin^2 w)` (elementary), with
`w = 2 pi / m`, give

> **P_5(regular m-gon) = 1 − 5 (15 cos²w + 92 cos w + 76) / (36 m² sin²w),  w = 2π/m.**

It returns Valtr's `11/36` at `m = 3`, Valtr's `49/144` at `m = 4`, and
`1 − 305/(48π²)` as `m → ∞`. We did not find this closed form in the literature, but it is a
one-line corollary of Alikoski + identity (II), so the novelty is only as strong as (II)'s
(see `LITERATURE.md` §5).

| m | $P_4$ (Sylvester) | $E[A_3^2]$ | $P_5$ (this session) | $P_5$ decimal |
|---|---|---|---|---|
| 3 | $\frac{2}{3}$ | $\frac{1}{72}$ | $\frac{11}{36}$ | 0.305555555555556 |
| 4 | $\frac{25}{36}$ | $\frac{1}{96}$ | $\frac{49}{144}$ | 0.340277777777778 |
| 5 | $\frac{4}{5} - \frac{2 \sqrt{5}}{45}$ | $\frac{\sqrt{5}}{1500} + \frac{1}{120}$ | $\frac{7}{12} - \frac{47 \sqrt{5}}{450}$ | 0.349788455683355 |
| 6 | $\frac{683}{972}$ | $\frac{25}{2592}$ | $\frac{1373}{3888}$ | 0.353137860082305 |
| 7 | 0.703522828397 | 0.009574699002 | 0.354554061016 | 0.354554061016242 |
| 8 | $\frac{479}{576} - \frac{13 \sqrt{2}}{144}$ | $\frac{\left(\sqrt{2} + 4\right)^{2}}{3072}$ | $\frac{1469}{2304} - \frac{115 \sqrt{2}}{576}$ | 0.355235139456761 |
| 9 | $\frac{-105 - 2952 \sin^{4}{\left(\frac{\pi}{9} \right)} + 3056 \sin^{2}{\left(\frac{\pi}{9} \right)}}{729 \sin^{2}{\left(\frac{2 \pi}{9} \right)}}$ | 0.009525484036 | 0.355595330567 | 0.355595330566635 |
| 10 | $\frac{751}{900} - \frac{131 \sqrt{5}}{2250}$ | $\frac{11 \sqrt{5}}{6000} + \frac{13}{2400}$ | $\frac{461}{720} - \frac{229 \sqrt{5}}{1800}$ | 0.355800240640305 |
| 11 | 0.704327295391 | 0.009510551536 | 0.355923753836 | 0.355923753835515 |
| 12 | $\frac{1093}{1296} - \frac{13 \sqrt{3}}{162}$ | $\frac{\left(\sqrt{3} + 4\right)^{2}}{3456}$ | $\frac{3439}{5184} - \frac{115 \sqrt{3}}{648}$ | 0.356001785693795 |
| $\infty$ (disk) | $1-\frac{35}{12\pi^2}$ | $\frac{3}{32\pi^2}$ | $1-\frac{305}{48\pi^2}$ | 0.356188312272645 |

`P_5` is strictly increasing in `m`, from the triangle to the disk — consistent with, and a
sharp quantitative form of, the extremal theorem restricted to regular polygons.

## 3. Numerical verification of the theorem over 561 bodies

`src/n5_landscape.py`, `results/n5_landscape.json`. Families, all evaluated **exactly**
(route P on the polygon, exact-`Fraction` covariance, affine-whitened first):

* regular `m`-gons, `m = 3..40`
* triangles with all three corners truncated, 25 values of the cut parameter
* circular caps (chord-cut segments), 25 opening angles from `0.05` to `2π`
* stadiums, 21 aspect ratios from `1e-3` to `8`
* Minkowski interpolations `(1−t)·triangle + t·disk`, 21 values of `t`
* regular `m`-gons with one vertex pushed out by factors `1.05 … 6`
* **400 random convex polygons** (hulls of Gaussian point clouds, 3–12 vertices)
* half-disk

**Result: 0 violations of the window** at tolerance `1e-12`.

| | F | P_5 | body |
|---|---|---|---|
| closest to the disk bound | 0.064381168788 | 0.356188312121 | 400-gon approximation of the disk (excess `1.5e-11` is the polygonal-approximation error) |
| closest from inside, genuine polygon | 0.064381318806 | 0.356186811936 | regular 40-gon |
| closest to the triangle bound | 0.069441850947 | 0.305581490531 | triangle with corners truncated at `t = 0.4990` |

Self-test run **before** the scan: `F` evaluated on 200 random affine images of a triangle
and 200 of a parallelogram (condition numbers up to `~1e6`) reproduces `5/72` and
`(1−49/144)/10` to `5.3e-16` and `7.2e-16`. This is the regression test that caught the bug
described in §5.

## 4. Direct optimisation over convex polygons

Nelder–Mead over the `2m` free coordinates, `F` evaluated on the convex hull, 6 restarts,
so the search ranges over all convex polygons with at most `m` vertices:

| search | F found | interpretation |
|---|---|---|
| max over ≤3,4,5,6,8,10-gons | 0.069444444444 in every case | the triangle bound `5/72`, attained — never exceeded |
| min over ≤4-gons | 0.065972222222 | `= (1 − 49/144)/10`: the parallelogram, exactly |
| min over ≤5-gons | 0.065021154432 | the regular pentagon, exactly |
| min over ≤6-gons | 0.064686213992 | the regular hexagon, exactly |
| min over ≤8-gons | 0.064476486054 | the regular octagon, exactly |
| min over ≤10-gons | 0.064440631395 | **not** the regular decagon (`0.064419975899`) — Nelder–Mead under-converged in 20 dimensions; reported as a limitation, not a result |

The maximiser is a triangle for every `m` and the minimisers are the regular `m`-gons (up to
the `m = 10` convergence failure) — exactly the shape the theorem predicts.

## 5. A bug this search found (and what it teaches)

The first run reported `max over ≤3-gons: F = 0.076256 > 5/72`. Since every triangle is an
affine image of every other and `F` is affine invariant, `F` is *constant* on triangles, so
this was a bug by construction. Nelder–Mead had found an almost-collinear triangle (area
`2.4e-5`, edges `4.1, 1.4, 2.7`) where the float64 covariance `E[xx^T] − mu mu^T` cancelled
catastrophically: `E[A_3^2]` came out `0.00708` instead of `1/72`. Route P itself was
healthy there. Fixed by computing the covariance in exact `Fraction` arithmetic and by
putting every polygon in isotropic position first; the offending triangle then returns
`5/72` to `1e-16`.

**Lesson:** an affine-invariant quantity is its own regression test. Evaluate it on random
affine images with large condition numbers and require invariance to machine precision — and
run that test *before* any search, because optimisers hunt for numerical noise.

## 6. On a proof route (asked for by TASK.md B3b)

Superseded by Marckert–Rahmani, but recorded for completeness. The reduction
`P_5 = 1 − 10 E[A_3(1−A_3)]` says the n = 5 statement is *equivalent* to a sharp two-sided
bound on `F = E[A_3] − (3/2) det Sigma/|K|^2`, i.e. to a single inequality between
Sylvester's functional and the planar isotropic constant. We did not find a
symmetrisation argument that settles it (`F` is not obviously monotone under Steiner
symmetrisation — `E[A_3]` decreases but `det Sigma/|K|^2` also decreases, so the difference
is not controlled), and we did not pursue it further once the theorem was found.


**Superseded 2026-08-19 06:30.** A symmetrisation argument DOES settle it: with vertical
sections the two moments change in the exact fibrewise ratio `Delta E[A_3^2]/Delta E[A_3] = b_2/V < 1`,
so `F` is monotone under Steiner symmetrisation and shaking after all. Proof (external, verified
here to 9 digits on the half-disk/disk gap): [`N5_PROOF.md`](N5_PROOF.md).
