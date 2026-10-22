"""Incremental Bowyer-Watson Delaunay triangulation with ghost triangles and exact predicates.

A "ghost" vertex (index -1) stands for the point at infinity.  Every hull edge a->b has a ghost
triangle (b, a, -1) on its outside, so the whole triangulation is a closed surface and the cavity
search needs no special cases for the convex hull.  Points are inserted in Hilbert curve order and
located by a straight line walk starting from the previous insertion.
"""
import math
from fractions import Fraction

from .predicates import incircle, orient2d

G = -1


class DegenerateInput(ValueError):
    """Raised when fewer than three distinct, non collinear points are supplied."""


def hilbert_key(x, y, order=16):
    n = 1 << order
    d = 0
    s = n >> 1
    while s:
        rx = 1 if x & s else 0
        ry = 1 if y & s else 0
        d += s * s * ((3 * rx) ^ ry)
        if ry == 0:
            if rx == 1:
                x, y = n - 1 - x, n - 1 - y
            x, y = y, x
        s >>= 1
    return d


class Delaunay:
    def __init__(self, points):
        pts, self.input_to_vertex = [], []
        seen = {}
        for p in points:
            key = (p[0], p[1])
            if key not in seen:
                seen[key] = len(pts)
                pts.append(key)
            self.input_to_vertex.append(seen[key])
        self.points = pts
        self.tris = {}      # id -> (a, b, c)
        self.edge = {}      # (u, v) -> id of the triangle that has u->v as a directed edge
        self._next = 0
        self._hint = None
        self._build()

    # ------------------------------------------------------------------ construction
    def _add(self, a, b, c):
        t = self._next
        self._next += 1
        self.tris[t] = (a, b, c)
        self.edge[(a, b)] = t
        self.edge[(b, c)] = t
        self.edge[(c, a)] = t
        if G not in (a, b, c):
            self._hint = t
        return t

    def _remove(self, t):
        a, b, c = self.tris.pop(t)
        del self.edge[(a, b)], self.edge[(b, c)], self.edge[(c, a)]

    def _build(self):
        pts = self.points
        n = len(pts)
        if n < 3:
            raise DegenerateInput("need at least 3 distinct points")
        i0, i1, i2 = 0, 1, None
        for k in range(2, n):
            if orient2d(pts[i0], pts[i1], pts[k]) != 0:
                i2 = k
                break
        if i2 is None:
            raise DegenerateInput("all points are collinear")
        if orient2d(pts[i0], pts[i1], pts[i2]) < 0:
            i1, i2 = i2, i1
        a, b, c = i0, i1, i2
        self._add(a, b, c)
        self._add(b, a, G)
        self._add(c, b, G)
        self._add(a, c, G)
        rest = [i for i in range(n) if i not in (a, b, c)]
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        x0, y0 = min(xs), min(ys)
        sx = (max(xs) - x0) or 1.0
        sy = (max(ys) - y0) or 1.0
        top = (1 << 16) - 1
        rest.sort(key=lambda i: hilbert_key(int((pts[i][0] - x0) / sx * top), int((pts[i][1] - y0) / sy * top)))
        for i in rest:
            self._insert(i)

    def _conflict(self, t, p):
        """Is p strictly inside the (generalised) circumcircle of triangle t?"""
        tri = self.tris[t]
        pts = self.points
        if G in tri:
            i = tri.index(G)
            u, v = tri[(i + 1) % 3], tri[(i + 2) % 3]   # ghost (u, v, G); hull edge is v->u
            pu, pv = pts[u], pts[v]
            o = orient2d(pv, pu, p)
            if o != 0:
                return o < 0
            # collinear with the hull edge: conflict only when strictly inside the segment
            return _strictly_between(pv, pu, p)
        return incircle(pts[tri[0]], pts[tri[1]], pts[tri[2]], p) > 0

    def _locate(self, p):
        """Return the id of some triangle in conflict with p."""
        pts = self.points
        t = self._hint
        for _ in range(4 * len(self.tris) + 16):
            tri = self.tris[t]
            if G in tri:
                if self._conflict(t, p):
                    return t
                break
            for k in range(3):
                u, v = tri[k], tri[(k + 1) % 3]
                if orient2d(pts[u], pts[v], p) < 0:
                    t = self.edge[(v, u)]
                    break
            else:
                return t
        for t in self.tris:  # fallback; essentially never taken on valid input
            if self._conflict(t, p):
                return t
        raise AssertionError("no conflicting triangle found")

    def _insert(self, idx):
        p = self.points[idx]
        seed = self._locate(p)
        cavity = {seed}
        stack = [seed]
        rejected = set()
        while stack:
            t = stack.pop()
            a, b, c = self.tris[t]
            for u, v in ((a, b), (b, c), (c, a)):
                nb = self.edge[(v, u)]
                if nb in cavity or nb in rejected:
                    continue
                if self._conflict(nb, p):
                    cavity.add(nb)
                    stack.append(nb)
                else:
                    rejected.add(nb)
        boundary = []
        for t in cavity:
            a, b, c = self.tris[t]
            for u, v in ((a, b), (b, c), (c, a)):
                if self.edge[(v, u)] not in cavity:
                    boundary.append((u, v))
        for t in cavity:
            self._remove(t)
        for u, v in boundary:
            self._add(u, v, idx)

    # ------------------------------------------------------------------ queries
    @property
    def triangles(self):
        """Real triangles as CCW vertex index triples, sorted for determinism."""
        return sorted(t for t in self.tris.values() if G not in t)

    def hull(self):
        """Convex hull vertex indices in counter clockwise order (collinear hull points included)."""
        nxt = {}
        for tri in self.tris.values():
            if G in tri:
                i = tri.index(G)
                u, v = tri[(i + 1) % 3], tri[(i + 2) % 3]
                nxt[v] = u
        start = min(nxt)
        out, cur = [start], nxt[start]
        while cur != start:
            out.append(cur)
            cur = nxt[cur]
        return out

    def neighbors(self):
        """Adjacency sets of the Delaunay graph (real vertices only)."""
        adj = {i: set() for i in range(len(self.points))}
        for a, b, c in self.triangles:
            adj[a].update((b, c))
            adj[b].update((a, c))
            adj[c].update((a, b))
        return adj

    def nearest(self, q, start=0, adj=None):
        """Exact nearest site by greedy routing on the Delaunay graph (correct for Delaunay graphs)."""
        adj = adj or self.neighbors()
        pts = self.points
        d2 = lambda i: (pts[i][0] - q[0]) ** 2 + (pts[i][1] - q[1]) ** 2
        cur, best = start, d2(start)
        improved = True
        while improved:
            improved = False
            for j in adj[cur]:
                dj = d2(j)
                if dj < best:
                    cur, best, improved = j, dj, True
        return cur

    def circumcenter(self, tri):
        (ax, ay), (bx, by), (cx, cy) = (self.points[i] for i in tri)
        ax, ay, bx, by, cx, cy = (Fraction(v) for v in (ax, ay, bx, by, cx, cy))
        d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
        ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d
        uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d
        return float(ux), float(uy)

    # ------------------------------------------------------------------ verification
    def validate(self, brute_force=True):
        """Check structural invariants; with brute_force also the empty circle property (O(n * t))."""
        n, tris = len(self.points), self.triangles
        h = len(self.hull())
        assert len(tris) == 2 * n - 2 - h, "Euler: triangle count"
        for t in tris:
            assert orient2d(*(self.points[i] for i in t)) > 0, "triangle not CCW"
        for tri in self.tris.values():
            a, b, c = tri
            for u, v in ((a, b), (b, c), (c, a)):
                assert (v, u) in self.edge, "missing twin edge"
        if brute_force:
            for a, b, c in tris:
                pa, pb, pc = self.points[a], self.points[b], self.points[c]
                for i, p in enumerate(self.points):
                    if i not in (a, b, c):
                        assert incircle(pa, pb, pc, p) <= 0, "empty circle violated"
        return True


def _strictly_between(a, b, p):
    """p is collinear with a,b; is it strictly inside the open segment ab?"""
    ax, ay, bx, by, px, py = (Fraction(v) for v in (a[0], a[1], b[0], b[1], p[0], p[1]))
    return (px - ax) * (bx - ax) + (py - ay) * (by - ay) > 0 and (px - bx) * (ax - bx) + (py - by) * (ay - by) > 0
