import argparse
import random
import sys
import time

from . import Delaunay, lloyd, voronoi_cells
from .svg import render


def main(argv=None):
    ap = argparse.ArgumentParser(prog="delaunay_voronoi", description="Delaunay triangulation and Voronoi diagram demo")
    ap.add_argument("-n", type=int, default=200, help="number of random points (ignored with --input)")
    ap.add_argument("--input", help="text file with one 'x y' pair per line")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--lloyd", type=int, default=0, help="Lloyd relaxation iterations")
    ap.add_argument("--svg", help="write an SVG picture to this path")
    ap.add_argument("--verify", action="store_true", help="run the brute force empty circle check")
    a = ap.parse_args(argv)
    if a.input:
        with open(a.input) as fh:
            pts = [tuple(map(float, ln.split()[:2])) for ln in fh if ln.strip()]
    else:
        rnd = random.Random(a.seed)
        pts = [(rnd.random(), rnd.random()) for _ in range(a.n)]
    bbox = (min(p[0] for p in pts), min(p[1] for p in pts), max(p[0] for p in pts), max(p[1] for p in pts))
    if a.lloyd:
        pts = lloyd(pts, bbox, a.lloyd)
    t = time.perf_counter()
    dt = Delaunay(pts)
    el = time.perf_counter() - t
    print("points=%d triangles=%d hull=%d time=%.3fs" % (len(dt.points), len(dt.triangles), len(dt.hull()), el))
    if a.verify:
        dt.validate()
        print("verified: Euler counts, orientation, empty circle property")
    if a.svg:
        with open(a.svg, "w") as fh:
            fh.write(render(dt, voronoi_cells(dt, bbox)))
        print("wrote", a.svg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
