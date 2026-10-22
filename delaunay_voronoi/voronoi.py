"""Bounded Voronoi diagram from the Delaunay graph by half plane clipping, and Lloyd relaxation."""
from .delaunay import Delaunay


def _clip(poly, nx, ny, c):
    """Keep the part of convex polygon poly with nx*x + ny*y <= c (Sutherland-Hodgman)."""
    out = []
    m = len(poly)
    for i in range(m):
        p, q = poly[i], poly[(i + 1) % m]
        fp, fq = nx * p[0] + ny * p[1] - c, nx * q[0] + ny * q[1] - c
        if fp <= 0:
            out.append(p)
        if (fp < 0 < fq) or (fq < 0 < fp):
            t = fp / (fp - fq)
            out.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    return out


def polygon_area(poly):
    return 0.5 * sum(poly[i][0] * poly[(i + 1) % len(poly)][1] - poly[(i + 1) % len(poly)][0] * poly[i][1]
                     for i in range(len(poly)))


def centroid(poly):
    a = polygon_area(poly)
    if abs(a) < 1e-300:
        return (sum(p[0] for p in poly) / len(poly), sum(p[1] for p in poly) / len(poly))
    cx = cy = 0.0
    for i in range(len(poly)):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % len(poly)]
        w = x0 * y1 - x1 * y0
        cx += (x0 + x1) * w
        cy += (y0 + y1) * w
    return cx / (6 * a), cy / (6 * a)


def voronoi_cells(dt: Delaunay, bbox):
    """Return {vertex: polygon} clipped to bbox = (xmin, ymin, xmax, ymax).

    Every Voronoi neighbour of a site is a Delaunay neighbour, so clipping the box by the
    bisectors with the Delaunay neighbours yields the exact cell."""
    x0, y0, x1, y1 = bbox
    adj = dt.neighbors()
    cells = {}
    for i, p in enumerate(dt.points):
        poly = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
        for j in adj[i]:
            q = dt.points[j]
            nx, ny = q[0] - p[0], q[1] - p[1]
            c = nx * (p[0] + q[0]) / 2 + ny * (p[1] + q[1]) / 2
            poly = _clip(poly, nx, ny, c)
            if not poly:
                break
        cells[i] = poly
    return cells


def lloyd(points, bbox, iterations=10):
    """Lloyd relaxation: move every site to the centroid of its Voronoi cell."""
    pts = [tuple(p) for p in points]
    for _ in range(iterations):
        dt = Delaunay(pts)
        cells = voronoi_cells(dt, bbox)
        pts = [centroid(cells[i]) if cells[i] else dt.points[i] for i in range(len(dt.points))]
    return pts
