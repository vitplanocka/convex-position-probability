# Literature-novelty check prompt (convex-position campaign)

Paste the block below into a research-capable model/session with web access.

---

You are a mathematician doing a rigorous prior-art / novelty check for a geometric-probability
project. I need to know, for each result below, whether it is ALREADY IN THE PUBLISHED LITERATURE
(with an exact citation and quote) or appears to be absent. Be skeptical and thorough; do not
fabricate citations, values, or page numbers. For every claim attach a source URL you actually
consulted, and label it [verified-from-source] (you read the statement in the primary source),
[secondary] (a reliable secondary source), or [not-found] (searched, absent). A "not found" with a
list of what you searched is a valuable answer.

## Setting and notation

K = a convex body in the plane, area V. Points are i.i.d. uniform in K. Define the DIMENSIONLESS
random-triangle area A_3 = area(triangle of 3 uniform points) / V, and more generally
A_k = area(convex hull of k uniform points) / V. I want to know the novelty status of these
EXACT results:

### Group 1 -- moments of the random-triangle area (the main question)
- E[A_3^3] (third moment) for:  triangle = 31/9000,  square = 137/72000,  disk = 1001/(6400 pi^4)

  Normalization anchors to confirm you have the right object: E[A_3] = 1/12 (triangle),
  11/144 (square), 35/(48 pi^2) (disk);  E[A_3^2] = 1/72 (triangle), 1/96 (square),
  3/(32 pi^2) (disk).  Question: is E[A_3^3] (for the triangle, square, or disk) published
  anywhere?  Also collect E[A_3^k] for k = 3, 4 if any source tabulates them.

### Group 2 -- second moment of the 4-point convex-hull area
- E[A_4^2] (second moment of the area of the convex hull of 4 uniform points):
  triangle = 181/4500, square = 859/27000, disk = (2400 pi^2 + 31031)/(19200 pi^4).
- Its "convex quadrilateral" part E[A_4^2 & the 4 points are in convex position]:
  triangle = 119/4500, square = 1307/54000, disk = (2400 pi^2 + 19019)/(19200 pi^4).
  Equivalently, the VARIANCE of the 4-point hull area, or Var of the convex-hull area of n points
  for small n, in the triangle/square/disk.  Is any of this published?

### Group 3 -- identities
- "E[area of convex hull of 4 uniform points] = 2 * E[area of a random triangle]"
  (i.e. E[A_4] = 2 E[A_3]) for EVERY planar distribution.  Is this stated anywhere as a general
  identity?  (Its Gaussian instance is a known curiosity; I want the general convex-body / general
  law statement.)
- "P(5 points in convex position) = 1 - 10 E[A_3 (1 - A_3)] = (5/2) P_4 - 3/2 + 15 det(Sigma)/V^2"
  where P_4 = Sylvester four-point probability and Sigma = covariance of the uniform law on K.
  Is this exact five-point formula (in terms of the first two triangle-area moments / a covariance)
  published?  Note: Marckert & Rahmani, "Around Sylvester's question in the plane" (Mathematika 67
  (2021) 860-884; arXiv:1511.03658) is the closest paper -- I have read it; its "new formula" is a
  different object (a comb recursion), and this identity is NOT in it.  Check other sources.

### Group 4 -- convex-position probabilities for regular polygons at n = 5, 6
- P_5(regular m-gon) = 1 - 5(15 cos^2 w + 92 cos w + 76)/(36 m^2 sin^2 w), w = 2 pi/m
  (gives Valtr's 11/36 at m=3, 49/144 at m=4, and 1 - 305/(48 pi^2) as m -> infinity).
- P_6 for regular polygons:  pentagon = 8941/22500 - 1349 sqrt(5)/11250 = 0.129248...,
  hexagon = 461299/3499200 = 0.131830..., octagon = 30103/61440 - 116141 sqrt(2)/460800.
  Valtr (1995/96) gave exact P_n only for the triangle and parallelogram; Marckert (2017) did the
  disk.  Has anyone published exact P_5 or P_6 for the regular pentagon / hexagon / general m-gon?

## Sources to check specifically (find, read, quote)
1. W. J. Reed, "Random points in a simplex", Pacific J. Math. 54 (1974) 183-198 -- moments of the
   volume of a random simplex in a simplex; for the 2-simplex extract E[A_3^k] and evaluate k=3,4.
   Quote the general moment formula.
2. J. Philip, "The area of a random triangle in a square" / "... in a disk" / "... in a triangle"
   (KTH reports, ~2005-2010; author page math.kth.se/~johanph) -- these tabulate area moments.
3. C. Buchta: "Zufallspolygone in konvexen Vielecken", J. reine angew. Math. 347 (1984);
   "An identity relating moments of functionals of convex hulls", Discrete Comput. Geom. 33 (2005)
   125-142; "Exact formulae for variances of functionals of convex hulls" (if it exists);
   "On the number of vertices of the convex hull of random points in a square and a triangle"
   (Anz. Oesterreich. Akad. Wiss. 2009/10). Extract any Var(hull area) or E[A_4^2] / E[A_3^3].
4. V. S. Alagar, "On the distribution of a random triangle", J. Appl. Probab. 14 (1977) 284-297.
5. A. M. Mathai, "An Introduction to Geometrical Probability" (1999); Kendall & Moran,
   "Geometrical Probability" (1963); H. Solomon, "Geometric Probability" (1978) -- moment tables.
6. Buchta & Reitzner on random polytope hull areas/variances; Finch's "Mathematical Constants"
   notes on geometric probability (he compiles such values).
7. OEIS: search the exact rationals (31/9000, 137/72000, 181/4500, 859/27000, 119/4500,
   1307/54000) and their decimals.

## Deliverable
A structured report: one section per Group above, each item marked KNOWN (with citation + exact
quote) / NOT FOUND (with the searches tried) / PARTIAL; a consolidated reference list with URLs;
and a one-paragraph bottom line: which of these results, if any, appear genuinely absent from the
literature, and which are already published (and where). Do not overclaim novelty -- "not found"
is the strongest allowed verdict.
