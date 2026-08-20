# Convex Position Probability

**How likely are _n_ random points to be in convex position?**

Draw _n_ i.i.d. uniform points from a convex body _K_ in the plane. Let
_P_K(n)_ be the probability that all _n_ are vertices of their convex hull
(i.e. no point falls inside the hull of the others). This is affine-invariant,
so "square" means every parallelogram and "disk" every ellipse.

This repository collects exact closed forms, an exact computational machine
that produces them, and a large numerical attack on the open extremal
conjecture — for _n_ = 5 and _n_ = 6 in the plane, and their _R^d_ analogues.
Every value is cross-checked by at least two independent methods.

## Headline results

Exact closed forms found here that a 2026 prior-art sweep did **not** locate in
the literature (see [`LITERATURE.md`](LITERATURE.md) for the per-result novelty
verdicts — "not found" means "not located in accessible prior art," not
"provably absent"):

> **Regular polygons, _n_ = 5.** For the regular _m_-gon, with _w_ = 2π/_m_,
> ```
> P_5(m-gon) = 1 - 5(15 cos^2 w + 92 cos w + 76) / (36 m^2 sin^2 w)
> ```
> One formula that reproduces Valtr's 11/36 (triangle, _m_=3) and 49/144
> (square, _m_=4), and tends to 1 - 305/(48π^2) as _m_ → ∞ (the disk).

> **Regular polygons, _n_ = 6.**
> ```
> P_6(pentagon) = 8941/22500 - 1349*sqrt(5)/11250   = 0.129248...
> P_6(hexagon)  = 461299/3499200                     = 0.131830...
> P_6(octagon)  = 30103/61440 - 116141*sqrt(2)/460800 = 0.133516...
> ```

> **The _n_ = 6 triangle-area / hull-area moments** for the triangle, square
> and disk — including the disk higher moments E[A_3^3] = 1001/(6400π^4) and
> the 4-point hull-area second moment E[A_4^2] — assembled into exact
> P_6(triangle) = 91/900, P_6(square) = 49/400, and
> P_6(disk) = 1 - (146400π^2 - 473473)/(11520π^4) = 0.134309...

> **A covariance form of the five-point identity.** For any convex body,
> ```
> P_5 = 1 - 10 E[A_3(1 - A_3)] = (5/2) P_4 - 3/2 + 15 det(Sigma)/|K|^2
> ```
> where P_4 is the Sylvester four-point probability and Sigma the covariance
> of the uniform law on _K_.

> **An _R^d_ moment identity in one line.** For _d_+2 points in general
> position in _R^d_, a Radon-partition argument gives, pointwise,
> vol conv(all) = ½ Σ_i vol conv(all − i), hence
> ```
> E[A_{d+2}] = ((d+2)/2) * E[A_{d+1}]   in every dimension
> ```
> (the _d_=2 case is the classical E[A_4] = 2 E[A_3]). This one is flagged as
> **folklore, not claimed new** — but the proof is clean and dimension-free.

## What is _not_ new (recorded honestly)

The triangle and square third moments E[A_3^3] = 31/9000 and 137/72000 are
**already published** (Reed 1974; Beck, arXiv:2412.07952) and are _not_ claimed
here. The E[A_4] = 2 E[A_3] identity and its _R^d_ generalization above are
distribution-independent Efron/Buchta-type folklore. See
[`LITERATURE.md`](LITERATURE.md).

## The open problem, and what we did about it

Blaschke (1917) proved the triangle minimizes and the ellipse maximizes the
Sylvester four-point probability; Marckert–Rahmani (2021) settled _n_ = 5. The
extremal conjecture for **_n_ ≥ 6 is open**: is _P_K(n)_ still minimized by the
triangle and maximized by the ellipse?

This repository does **not** prove it. It supplies evidence: an exact-machine +
Monte-Carlo scan over **568 planar bodies** (_n_ = 6 and 7) and **240 bodies in
_R^3_** (_n_ = 5) found **zero counterexamples**, with the simplex a measured
strict local minimum and the ball a strict local maximum. That is consistent
with the conjecture, and no more. Full account in [`VERDICT.md`](VERDICT.md).

## Contents

| Path | What's in it |
|---|---|
| [`THEOREMS.md`](THEOREMS.md) | Statements, proofs, and the two independent machine checks behind every exact value |
| [`LITERATURE.md`](LITERATURE.md) | Prior-art review and the per-result novelty verdicts |
| [`VERDICT.md`](VERDICT.md) | The extremal-conjecture scan: method, bodies tested, and the honest "consistent, not proved" conclusion |
| [`PROGRESS.md`](PROGRESS.md) | Chronological investigation log |
| [`N6_INGREDIENTS.md`](N6_INGREDIENTS.md) | The _n_ = 6 moment ingredient table |
| [`docs/`](docs/) | Extended write-ups: the _n_ = 5 proof (`N5_PROOF.md`), the _n_ = 6 landscape and extremal analysis, the _R^d_ analogues (`RD_ANALOGUES.md`), and the literature-check materials |
| [`src/`](src/) | Python: exact integrators (width-function and two-chord Blaschke–Petkantschin routes), samplers/hull testers, moment and scan code |
| [`results/`](results/) | JSON/log outputs behind the tables and scans |
| [`explainer/`](explainer/) | A self-contained HTML explainer (`in-convex-position.html`) — open in a browser, no build step |

## Quick start

```
pip install -r requirements.txt

python src/anchors.py              # exact anchors: Valtr formulas, Sylvester disk value
python src/n6_bp_polygon.py        # width-function route to n=6 polygon moments
python src/n6_twochord_polygon.py  # independent two-chord route (cross-check)
python src/rd_identity.py          # R^d moment identity, Monte-Carlo confirmation
```

The two `n6_*` routes are deliberately independent implementations; agreement
between them (and with Monte Carlo, and with every Valtr/Marckert closed-form
anchor) is what earns each exact value its place in `THEOREMS.md`.

## Status

The exact closed forms and the two-route machine are done and cross-verified.
The extremal conjecture for _n_ ≥ 6 remains open; this repository adds strong,
consistent numerical evidence but no proof. The one residual novelty risk is a
single unretrievable prior source (Philip's KTH report `area12.pdf`, currently
404), noted in `LITERATURE.md`.

## License

MIT — see [`LICENSE`](LICENSE).
