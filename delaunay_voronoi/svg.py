"""Tiny SVG renderer for a triangulation and its Voronoi cells."""


def render(dt, cells=None, size=700, pad=20):
    xs = [p[0] for p in dt.points]
    ys = [p[1] for p in dt.points]
    x0, y0 = min(xs), min(ys)
    span = max(max(xs) - x0, max(ys) - y0) or 1.0
    k = (size - 2 * pad) / span
    f = lambda p: "%.2f,%.2f" % (pad + (p[0] - x0) * k, size - pad - (p[1] - y0) * k)
    out = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">' % (size, size, size, size),
           '<rect width="100%" height="100%" fill="white"/>']
    if cells:
        for poly in cells.values():
            if poly:
                out.append('<polygon points="%s" fill="none" stroke="#d9534f" stroke-width="1"/>' % " ".join(map(f, poly)))
    for a, b, c in dt.triangles:
        out.append('<polygon points="%s" fill="none" stroke="#4477aa" stroke-width="0.8"/>' % " ".join(f(dt.points[i]) for i in (a, b, c)))
    for p in dt.points:
        out.append('<circle cx="%s" cy="%s" r="2" fill="black"/>' % tuple(f(p).split(",")))
    out.append("</svg>")
    return "\n".join(out)
