r"""Extra convex bodies for the identity tests (A3): long ellipse, half-disk, stadium,
plus their EXACT covariance determinants (hence exact E[A_3^2] = (3/2) det Sigma/|K|^2),
and polygonal approximations for route P.

Bodies (all convex):
  ellipse3   x^2/9 + y^2 <= 1        (affine image of the disk: every quantity must
                                      equal the disk's -- a null test of the pipeline)
  halfdisk   x^2 + y^2 <= 1, y >= 0
  stadium    convex hull of the two unit disks centred at (-1,0) and (1,0)
  cutsquare  unit square with one corner cut (a non-symmetric polygon, exact rational)
"""
import math
import numpy as np


# ------------------------------------------------------------------- samplers

def sample_ellipse3(rng, m, n):
    r = np.sqrt(rng.random((m, n)))
    t = rng.random((m, n)) * (2 * np.pi)
    return np.stack([3 * r * np.cos(t), r * np.sin(t)], axis=-1)


def sample_halfdisk(rng, m, n):
    r = np.sqrt(rng.random((m, n)))
    t = rng.random((m, n)) * np.pi
    return np.stack([r * np.cos(t), r * np.sin(t)], axis=-1)


def sample_stadium(rng, m, n, a=1.0):
    """Stadium = [-a,a]x[-1,1] rectangle plus unit half-disks at x = +-a. Area = 4a + pi."""
    area_rect = 4 * a
    area_tot = area_rect + math.pi
    u = rng.random((m, n))
    inrect = u < area_rect / area_tot
    x = np.empty((m, n))
    y = np.empty((m, n))
    # rectangle part
    x_r = rng.random((m, n)) * (2 * a) - a
    y_r = rng.random((m, n)) * 2 - 1
    # disk part (unit disk, then shifted left/right by a with sign of x)
    rr = np.sqrt(rng.random((m, n)))
    tt = rng.random((m, n)) * (2 * np.pi)
    xd = rr * np.cos(tt)
    yd = rr * np.sin(tt)
    xd = xd + np.where(xd >= 0, a, -a)
    x = np.where(inrect, x_r, xd)
    y = np.where(inrect, y_r, yd)
    return np.stack([x, y], axis=-1)


# -------------------------------------------------- exact covariance / E[A_3^2]

def exact_cov(body):
    """(det Sigma, |K|) exactly (as floats built from exact closed forms)."""
    pi = math.pi
    if body == "disk":
        return (1 / 4) ** 2, pi                      # Sigma = (1/4) I
    if body == "ellipse3":
        # affine image with matrix diag(3,1): det Sigma -> 9 det Sigma, |K| -> 3|K|
        return 9 * (1 / 4) ** 2, 3 * pi              # det Sigma/|K|^2 invariant
    if body == "halfdisk":
        A = pi / 2
        vx = 1 / 4
        vy = 1 / 4 - (4 / (3 * pi)) ** 2
        return vx * vy, A
    if body == "stadium":
        a = 1.0
        A = 4 * a + pi
        # E[y^2]: rectangle contributes int_{-a}^{a}int_{-1}^{1} y^2 = 2a*(2/3) = 4a/3
        #         two half-disks glue to one unit disk: int y^2 = pi/4
        Ey2 = (4 * a / 3 + pi / 4) / A
        # E[x^2]: rectangle int_{-a}^{a} x^2 dx * 2 = 4a^3/3.
        #   Right cap = {(a+u,v): u^2+v^2<=1, u>=0}: int (a+u)^2 = a^2 pi/2 + 2a*(2/3) + pi/8
        #   (int_{u>=0} u du dv = 2/3, int_{u>=0} u^2 = pi/8).  Left cap gives the same by
        #   x -> -x, so the two caps contribute a^2 pi + 8a/3 + pi/4.  (An earlier version
        #   dropped the 8a/3 term -- caught by the MC covariance check.)
        Ex2 = (4 * a ** 3 / 3 + a ** 2 * pi + 8 * a / 3 + pi / 4) / A
        return Ex2 * Ey2, A
    raise KeyError(body)


def exact_EA3sq(body):
    d, A = exact_cov(body)
    return 1.5 * d / A ** 2


# ------------------------------------------------- polygonal approximations

def poly_disk(m):
    return [(math.cos(2 * math.pi * k / m), math.sin(2 * math.pi * k / m)) for k in range(m)]


def poly_halfdisk(m):
    """m points on the semicircle (including both ends) -> convex polygon, CCW."""
    pts = [(math.cos(math.pi * k / (m - 1)), math.sin(math.pi * k / (m - 1))) for k in range(m)]
    # CCW order: start at (1,0), go over the top to (-1,0), then straight back
    return pts


def poly_stadium(m, a=1.0):
    """m points on each half-circle; CCW starting at (a,-1)."""
    right = [(a + math.cos(-math.pi / 2 + math.pi * k / (m - 1)),
              math.sin(-math.pi / 2 + math.pi * k / (m - 1))) for k in range(m)]
    left = [(-a + math.cos(math.pi / 2 + math.pi * k / (m - 1)),
             math.sin(math.pi / 2 + math.pi * k / (m - 1))) for k in range(m)]
    return right + left


def poly_ellipse3(m):
    return [(3 * math.cos(2 * math.pi * k / m), math.sin(2 * math.pi * k / m)) for k in range(m)]


POLY = {"disk": poly_disk, "halfdisk": poly_halfdisk, "stadium": poly_stadium,
        "ellipse3": poly_ellipse3}
SAMPLERS = {"ellipse3": sample_ellipse3, "halfdisk": sample_halfdisk, "stadium": sample_stadium}
AREAS = {"ellipse3": 3 * math.pi, "halfdisk": math.pi / 2, "stadium": 4 + math.pi}
