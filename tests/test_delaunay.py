import itertools
import math
import random
import unittest
from fractions import Fraction

from delaunay_voronoi import Delaunay, DegenerateInput, lloyd, orient2d, incircle, polygon_area, voronoi_cells
from delaunay_voronoi.predicates import incircle_exact, orient2d_exact, stats


def monotone_hull(pts):
    pts = sorted(set(pts))
    def half(seq):
        h = []
        for p in seq:
            while len(h) >= 2 and orient2d(h[-2], h[-1], p) <= 0:
                h.pop()
            h.append(p)
        return h[:-1]
    return half(pts) + half(reversed(pts))


class PredicateTests(unittest.TestCase):
    def test_orient_basic(self):
        self.assertEqual(orient2d((0, 0), (1, 0), (0, 1)), 1)
        self.assertEqual(orient2d((0, 0), (1, 0), (0, -1)), -1)
        self.assertEqual(orient2d((0, 0), (1, 1), (2, 2)), 0)

    def test_orient_nearly_collinear_matches_exact(self):
        rnd = random.Random(3)
        before = stats["orient_exact"]
        for _ in range(2000):
            x = rnd.random()
            a, b = (0.5 + x * 1e-9, 0.5), (12.0, 12.0)
            c = (24.0 + rnd.choice([0, 1e-14, -1e-14]), 24.0)
            self.assertEqual(orient2d(a, b, c), orient2d_exact(a, b, c))
        self.assertGreater(stats["orient_exact"], before)  # the exact path really ran

    def test_incircle_basic_and_cocircular(self):
        a, b, c = (0, 0), (1, 0), (0, 1)
        self.assertEqual(incircle(a, b, c, (1, 1)), 0)           # cocircular
        self.assertEqual(incircle(a, b, c, (0.5, 0.5)), 1)
        self.assertEqual(incircle(a, b, c, (2, 2)), -1)

    def test_incircle_random_matches_exact(self):
        rnd = random.Random(4)
        for _ in range(3000):
            p = [(rnd.randint(-5, 5) * 0.1, rnd.randint(-5, 5) * 0.1) for _ in range(4)]
            if orient2d(*p[:3]) > 0:
                self.assertEqual(incircle(*p), incircle_exact(*p))


class DelaunayTests(unittest.TestCase):
    def test_random_uniform(self):
        for seed in range(5):
            rnd = random.Random(seed)
            pts = [(rnd.random(), rnd.random()) for _ in range(150)]
            dt = Delaunay(pts)
            self.assertTrue(dt.validate())
            self.assertEqual(sorted(dt.points[i] for i in dt.hull()), sorted(monotone_hull(pts)))

    def test_grid_is_fully_degenerate(self):
        pts = [(x, y) for x in range(9) for y in range(7)]
        random.Random(1).shuffle(pts)
        dt = Delaunay(pts)
        self.assertTrue(dt.validate())
        self.assertEqual(len(dt.triangles), 2 * 8 * 6)

    def test_cocircular_points(self):
        pts = [(math.cos(2 * math.pi * k / 12), math.sin(2 * math.pi * k / 12)) for k in range(12)]
        pts.append((0.0, 0.0))
        self.assertTrue(Delaunay(pts).validate())

    def test_duplicates_are_merged(self):
        pts = [(0, 0), (1, 0), (0, 1), (1, 1), (0, 0), (1, 1)]
        dt = Delaunay(pts)
        self.assertEqual(len(dt.points), 4)
        self.assertEqual(dt.input_to_vertex[4], dt.input_to_vertex[0])
        self.assertTrue(dt.validate())

    def test_collinear_prefix_then_point(self):
        pts = [(0, 0), (1, 0), (2, 0), (3, 0), (1.5, 2)]
        dt = Delaunay(pts)
        self.assertTrue(dt.validate())
        self.assertEqual(len(dt.triangles), 3)

    def test_point_on_hull_edge(self):
        dt = Delaunay([(0, 0), (4, 0), (0, 4), (2, 0)])
        self.assertTrue(dt.validate())
        self.assertEqual(len(dt.triangles), 2)

    def test_degenerate_inputs(self):
        with self.assertRaises(DegenerateInput):
            Delaunay([(0, 0), (1, 1), (2, 2), (3, 3)])
        with self.assertRaises(DegenerateInput):
            Delaunay([(0, 0), (1, 1)])

    def test_clustered_and_tiny_scale(self):
        rnd = random.Random(9)
        pts = [(1e8 + rnd.random() * 1e-3, 1e8 + rnd.random() * 1e-3) for _ in range(80)]
        self.assertTrue(Delaunay(pts).validate())

    def test_nearest_matches_brute_force(self):
        rnd = random.Random(5)
        pts = [(rnd.random(), rnd.random()) for _ in range(200)]
        dt = Delaunay(pts)
        adj = dt.neighbors()
        for _ in range(100):
            q = (rnd.random() * 1.2 - 0.1, rnd.random() * 1.2 - 0.1)
            best = min(range(len(dt.points)), key=lambda i: (dt.points[i][0] - q[0]) ** 2 + (dt.points[i][1] - q[1]) ** 2)
            got = dt.nearest(q, start=rnd.randrange(len(dt.points)), adj=adj)
            self.assertAlmostEqual(math.dist(dt.points[got], q), math.dist(dt.points[best], q))

    def test_circumcenter_equidistant(self):
        dt = Delaunay([(0, 0), (4, 0), (0, 3), (4, 3)])
        for t in dt.triangles:
            c = dt.circumcenter(t)
            d = [math.dist(c, dt.points[i]) for i in t]
            self.assertAlmostEqual(d[0], d[1])
            self.assertAlmostEqual(d[1], d[2])


class VoronoiTests(unittest.TestCase):
    def test_cells_partition_the_box(self):
        rnd = random.Random(2)
        pts = [(rnd.random(), rnd.random()) for _ in range(120)]
        dt = Delaunay(pts)
        cells = voronoi_cells(dt, (-0.5, -0.5, 1.5, 1.5))
        self.assertAlmostEqual(sum(polygon_area(c) for c in cells.values()), 4.0, places=9)

    def test_cell_points_are_closest_to_own_site(self):
        rnd = random.Random(6)
        pts = [(rnd.random(), rnd.random()) for _ in range(60)]
        dt = Delaunay(pts)
        cells = voronoi_cells(dt, (0, 0, 1, 1))
        for i, poly in cells.items():
            cx = sum(p[0] for p in poly) / len(poly)
            cy = sum(p[1] for p in poly) / len(poly)
            near = min(range(len(dt.points)), key=lambda j: math.dist(dt.points[j], (cx, cy)))
            self.assertEqual(near, i)

    def test_lloyd_reduces_area_variance(self):
        rnd = random.Random(8)
        pts = [(rnd.random(), rnd.random()) for _ in range(60)]
        box = (0, 0, 1, 1)
        def var(p):
            areas = [polygon_area(c) for c in voronoi_cells(Delaunay(p), box).values()]
            m = sum(areas) / len(areas)
            return sum((a - m) ** 2 for a in areas)
        self.assertLess(var(lloyd(pts, box, 8)), var(pts))


if __name__ == "__main__":
    unittest.main()
