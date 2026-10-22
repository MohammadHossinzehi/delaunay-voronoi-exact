# delaunay-voronoi

Exact Delaunay triangulation and Voronoi diagrams in pure Python, with zero dependencies.

Delaunay code is notorious for failing on real data: grids, cocircular points and nearly collinear points make floating point predicates disagree with each other, and the algorithm then loops or produces garbage. This project makes the predicates exact, so those inputs just work, and it ships a test suite that attacks exactly those cases.

## Features

- **Robust predicates** (`predicates.py`): orient2d and incircle use a fast floating point evaluation with a rigorous error bound, and fall back to exact rational arithmetic (fractions.Fraction) only when the float result cannot be trusted. The returned sign is always the exact sign.
- **Bowyer-Watson Delaunay** (delaunay.py): incremental insertion with *ghost triangles* (a vertex at infinity), so the convex hull needs no special cases. Points are inserted in Hilbert curve order and located with a straight line walk from the previous insertion, which gives near linear behaviour on typical inputs.
- **Degeneracy handling**: duplicate points are merged, collinear prefixes and points on hull edges are handled, fully degenerate input raises DegenerateInput.
- **Queries**: convex hull, Delaunay graph, exact circumcenters, nearest site by greedy routing on the Delaunay graph.
- **Voronoi cells** (voronoi.py): bounded cells clipped to a box, built from the fact that Voronoi neighbours are Delaunay neighbours. Includes Lloyd relaxation.
- **Verification**: Delaunay.validate() checks Euler's formula (T = 2n - 2 - h), triangle orientation, edge twins and, by brute force, the empty circle property.
- **SVG output** for a quick picture.

## Run it

Requires Python 3.9+.

```bash
python -m delaunay_voronoi -n 3000 --verify
python -m delaunay_voronoi -n 150 --lloyd 10 --svg out.svg
python -m delaunay_voronoi --input points.txt --svg out.svg   # one "x y" pair per line
```

As a library:

```python
from delaunay_voronoi import Delaunay, voronoi_cells

dt = Delaunay([(0, 0), (4, 0), (0, 3), (4, 3), (2, 1)])
print(dt.triangles)      # CCW index triples
print(dt.hull())         # hull vertex indices, CCW
cells = voronoi_cells(dt, (-1, -1, 5, 4))
```

## Tests

```bash
python -m unittest discover -s tests
```

The suite covers: predicate signs against the exact evaluator on adversarial near degenerate input (and asserts the exact path really runs), uniform random sets checked against a monotone chain hull, fully degenerate integer grids, 12 cocircular points plus a centre, duplicates, collinear prefixes, points on hull edges, clustered points at 1e8 offset with 1e-3 spread, nearest neighbour against brute force, Voronoi cells partitioning the box exactly, every cell point being closest to its own site, and Lloyd relaxation lowering cell area variance.

## Design notes

- **Why ghost triangles?** Every hull edge a->b gets a ghost triangle (b, a, inf). A ghost is "in conflict" with a point when the point is strictly outside the hull edge, or strictly inside the segment. With that one rule the cavity search, the boundary extraction and the re-triangulation are identical for hull and interior insertions.
- **Why a float filter?** Exact rational arithmetic on every test is slow. The filter bounds (Shewchuk's) make the common case a handful of float operations.
- **Voronoi by clipping** is O(deg) per cell and trivially handles unbounded cells, at the cost of floats in the output geometry. The combinatorial structure (the triangulation) is exact.
- **Not included**: constrained Delaunay, deletion, and Shewchuk's multi stage adaptive expansions (the fallback here is simple exact rationals).
