"""Robust geometric predicates: a floating point filter with an exact rational fallback.

orient2d(a, b, c) > 0  <=> c is strictly left of the directed line a->b (a,b,c counter clockwise)
incircle(a, b, c, d) > 0 <=> d is strictly inside the circle through a,b,c (a,b,c counter clockwise)

The sign returned is always the exact sign for the given binary floating point inputs.
"""
from fractions import Fraction

_EPS = 2.0 ** -53
_ORIENT_BOUND = (3.0 + 16.0 * _EPS) * _EPS
_INCIRCLE_BOUND = (10.0 + 96.0 * _EPS) * _EPS

stats = {"orient_exact": 0, "incircle_exact": 0}


def _sign(x):
    return (x > 0) - (x < 0)


def orient2d(a, b, c):
    l = (a[0] - c[0]) * (b[1] - c[1])
    r = (a[1] - c[1]) * (b[0] - c[0])
    det = l - r
    if abs(det) > _ORIENT_BOUND * (abs(l) + abs(r)):
        return _sign(det)
    return orient2d_exact(a, b, c)


def orient2d_exact(a, b, c):
    stats["orient_exact"] += 1
    ax, ay, bx, by, cx, cy = (Fraction(v) for v in (a[0], a[1], b[0], b[1], c[0], c[1]))
    return _sign((ax - cx) * (by - cy) - (ay - cy) * (bx - cx))


def incircle(a, b, c, d):
    adx, ady = a[0] - d[0], a[1] - d[1]
    bdx, bdy = b[0] - d[0], b[1] - d[1]
    cdx, cdy = c[0] - d[0], c[1] - d[1]
    alift = adx * adx + ady * ady
    blift = bdx * bdx + bdy * bdy
    clift = cdx * cdx + cdy * cdy
    t1 = alift * (bdx * cdy - cdx * bdy)
    t2 = blift * (cdx * ady - adx * cdy)
    t3 = clift * (adx * bdy - bdx * ady)
    det = t1 + t2 + t3
    perm = (
        (abs(bdx * cdy) + abs(cdx * bdy)) * alift
        + (abs(cdx * ady) + abs(adx * cdy)) * blift
        + (abs(adx * bdy) + abs(bdx * ady)) * clift
    )
    if abs(det) > _INCIRCLE_BOUND * perm:
        return _sign(det)
    return incircle_exact(a, b, c, d)


def incircle_exact(a, b, c, d):
    stats["incircle_exact"] += 1
    dx, dy = Fraction(d[0]), Fraction(d[1])
    adx, ady = Fraction(a[0]) - dx, Fraction(a[1]) - dy
    bdx, bdy = Fraction(b[0]) - dx, Fraction(b[1]) - dy
    cdx, cdy = Fraction(c[0]) - dx, Fraction(c[1]) - dy
    det = (
        (adx * adx + ady * ady) * (bdx * cdy - cdx * bdy)
        + (bdx * bdx + bdy * bdy) * (cdx * ady - adx * cdy)
        + (cdx * cdx + cdy * cdy) * (adx * bdy - bdx * ady)
    )
    return _sign(det)
