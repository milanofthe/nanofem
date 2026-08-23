#!/usr/bin/env python3
"""Regenerate the report data files from the release binary.

Run from the repository root after `cargo build --release`:

    python3 report/data/gen.py

Writes cond.dat, patch.dat and the field grids and contours next to this
file.
No third party packages are used.
"""

import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
BIN = os.path.join(ROOT, "target", "release", "nanofem")
MESH = os.path.join(ROOT, "models", "patch.msh")
SUB = "mat sub eps 2.2 tand 0.001\npec pec\nabc open\nport 1 feed 0 0 1 200\n"


def run(deck, keep_stderr=False):
    path = os.path.join(HERE, "_tmp.nfm")
    with open(path, "w") as f:
        f.write(deck)
    p = subprocess.run([BIN, path], capture_output=True, text=True)
    os.remove(path)
    if p.returncode != 0:
        sys.exit("nanofem failed: " + p.stderr)
    return p.stderr if keep_stderr else p.stdout


def conditioning():
    """Pivot spread against frequency, with and without equilibration."""
    freqs = ["1e6", "3e6", "1e7", "3e7", "1e8", "3e8", "1e9", "2.4e9", "5e9", "1e10"]
    rows = []
    for f in freqs:
        err = run("mesh %s\n%ssweep lin %s %s 1\n" % (MESH, SUB, f, f), keep_stderr=True)
        line = [l for l in err.splitlines() if "pivot spread" in l][0]
        rows.append((f, line.split("spread ")[1].split(",")[0]))
    # the unequilibrated column is measured with equilibration disabled in the
    # solver and is kept here as recorded, since the current binary cannot
    # produce it
    plain = ["1.1e9", "1.2e8", "1.1e7", "1.2e6", "1.1e5", "1.2e4",
             "1.1e3", "8.4e2", "3.4e3", "5.8e5"]
    with open(os.path.join(HERE, "cond.dat"), "w") as f:
        f.write("f cond plain\n")
        for (fr, c), p in zip(rows, plain):
            f.write("%s %s %s\n" % (fr, c, p))


def sweep():
    """Reflection of the reference antenna over the resonance."""
    out = run("mesh %s\n%ssweep lin 2.1e9 2.8e9 36\n" % (MESH, SUB))
    with open(os.path.join(HERE, "patch.dat"), "w") as f:
        f.write("f db\n")
        for line in out.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            v = line.split()
            m = math.hypot(float(v[1]), float(v[2]))
            f.write("%.6e %.4f\n" % (float(v[0]), 20 * math.log10(max(m, 1e-12))))


def read_vtk(path):
    """Cell centroids and |E| of a legacy VTK written by nanofem."""
    ls = open(path).read().splitlines()
    i = next(k for k, l in enumerate(ls) if l.startswith("POINTS"))
    npt = int(ls[i].split()[1])
    pts = [tuple(map(float, ls[i + 1 + k].split())) for k in range(npt)]
    i = next(k for k, l in enumerate(ls) if l.startswith("CELLS"))
    ncell = int(ls[i].split()[1])
    cells = [tuple(map(int, ls[i + 1 + k].split()[1:])) for k in range(ncell)]
    i = next(k for k, l in enumerate(ls) if l.startswith("SCALARS Emag"))
    mag = [float(ls[i + 2 + k]) for k in range(ncell)]
    cen = []
    for c in cells:
        cen.append(tuple(sum(pts[v][d] for v in c) / 4.0 for d in range(3)))
    return cen, mag


def grid_idw(cen, mag, axis, lo, hi, u, v, nu, nv, bounds, radius, floor_db):
    """Smooth |E| in dB onto a regular grid, gaussian weighted.

    The field is piecewise constant per tetrahedron, so a raster of the raw
    cell values is blocky at any useful resolution. Weighting all cells of
    the slab by exp(-(d/r)^2) gives a smooth field whose support is the mesh
    spacing. The result is in dB relative to the peak of the cut, since the
    fringing field above a patch is two orders of magnitude below the field
    inside the substrate.
    """
    u0, u1, v0, v1 = bounds
    pts = [(c[u], c[v], m) for c, m in zip(cen, mag) if lo <= c[axis] <= hi]
    if not pts:
        sys.exit("empty slab")
    out = []
    for r in range(nv):
        y = v0 + (v1 - v0) * r / (nv - 1)
        row = []
        for k in range(nu):
            x = u0 + (u1 - u0) * k / (nu - 1)
            num = den = 0.0
            for cu, cv, m in pts:
                w = math.exp(-(((cu - x) ** 2 + (cv - y) ** 2) / (radius * radius)))
                num += w * m
                den += w
            row.append(num / den if den > 0 else 0.0)
        out.append(row)
    peak = max(max(r) for r in out) or 1.0
    lin = 10 ** (floor_db / 20.0)
    return [[max(20.0 * math.log10(max(m / peak, lin)), floor_db) for m in r] for r in out]


def contours(acc, bounds, nu, nv, scale, levels):
    """Marching squares on the normalized grid, as polyline segments."""
    u0, u1, v0, v1 = bounds
    segs = []
    for lv in levels:
        for r in range(nv - 1):
            for k in range(nu - 1):
                cs = [acc[r][k], acc[r][k + 1], acc[r + 1][k + 1], acc[r + 1][k]]
                xs = [u0 + (u1 - u0) * (k + dx) / (nu - 1) for dx in (0, 1, 1, 0)]
                ys = [v0 + (v1 - v0) * (r + dy) / (nv - 1) for dy in (0, 0, 1, 1)]
                hit = []
                for e in range(4):
                    a, b = e, (e + 1) % 4
                    if (cs[a] - lv) * (cs[b] - lv) < 0:
                        t = (lv - cs[a]) / (cs[b] - cs[a])
                        hit.append((xs[a] + t * (xs[b] - xs[a]), ys[a] + t * (ys[b] - ys[a])))
                if len(hit) == 2:
                    segs.append((hit[0], hit[1]))
    return [((a[0] * scale, a[1] * scale), (b[0] * scale, b[1] * scale)) for a, b in segs]


def write_grid(name, acc, bounds, nu, nv, scale, floor_db):
    u0, u1, v0, v1 = bounds
    with open(os.path.join(HERE, name), "w") as f:
        f.write("x y e\n")
        for r in range(nv):
            for k in range(nu):
                x = (u0 + (u1 - u0) * k / (nu - 1)) * scale
                y = (v0 + (v1 - v0) * r / (nv - 1)) * scale
                f.write("%.4f %.4f %.4f\n" % (x, y, (acc[r][k] - floor_db) / (-floor_db)))
            f.write("\n")


def write_contours(name, segs):
    with open(os.path.join(HERE, name), "w") as f:
        f.write("x y\n")
        for a, b in segs:
            f.write("%.4f %.4f\n%.4f %.4f\n\n" % (a[0], a[1], b[0], b[1]))


def mesh_cut(plane=0.045, axis=0, u=1, v=2, name="mesh_cut.dat"):
    """Intersection of the tetrahedral mesh with a plane, as polygons.

    Each tetrahedron crossing the plane is cut into a triangle or a
    quadrilateral. The polygon vertices are ordered by angle about their
    centroid, which is enough for a convex cross section.
    """
    ls = open(MESH).read().splitlines()
    i = ls.index("$Nodes")
    n = int(ls[i + 1])
    pos = {}
    for k in range(n):
        t = ls[i + 2 + k].split()
        pos[int(t[0])] = (float(t[1]), float(t[2]), float(t[3]))
    i = ls.index("$Elements")
    tets = []
    for k in range(int(ls[i + 1])):
        t = ls[i + 2 + k].split()
        if t[1] == "4":
            nt = int(t[2])
            tets.append([int(x) for x in t[3 + nt:]])
    ED = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    polys = []
    for tet in tets:
        pts = []
        for a, b in ED:
            pa, pb = pos[tet[a]], pos[tet[b]]
            da, db = pa[axis] - plane, pb[axis] - plane
            if da * db < 0:
                w = da / (da - db)
                pts.append((pa[u] + w * (pb[u] - pa[u]), pa[v] + w * (pb[v] - pa[v])))
        if len(pts) < 3:
            continue
        cx = sum(q[0] for q in pts) / len(pts)
        cy = sum(q[1] for q in pts) / len(pts)
        pts.sort(key=lambda q: math.atan2(q[1] - cy, q[0] - cx))
        polys.append(pts)
    with open(os.path.join(HERE, name), "w") as f:
        f.write("x y\n")
        for q in polys:
            for a in q + [q[0]]:
                f.write("%.4f %.4f\n" % (a[0] * 1000.0, a[1] * 1000.0))
            f.write("\n")
    return len(polys)


def fields():
    """Two cuts of |E| through the antenna at its resonance."""
    vtk = os.path.join(HERE, "_field.vtk")
    run("mesh %s\n%ssweep lin 2.45e9 2.45e9 1\nfield %s 2.45e9\n" % (MESH, SUB, vtk))
    cen, mag = read_vtk(vtk)
    os.remove(vtk)
    h = 0.00157
    lv = [0.35, 0.55, 0.75, 0.9]
    # inside the substrate, looking down on the patch
    bx, nu, nv = (0.012, 0.078, 0.022, 0.078), 90, 78
    acc = grid_idw(cen, mag, 2, 0.0, h, 0, 1, nu, nv, bx, 0.005, -20.0)
    write_grid("field_xy.dat", acc, bx, nu, nv, 1000.0, -20.0)
    nrm = [[(m + 20.0) / 20.0 for m in r] for r in acc]
    write_contours("cont_xy.dat", contours(nrm, bx, nu, nv, 1000.0, lv))
    # vertical cut along the resonant direction, through the middle of the
    # patch: the half wave in the substrate and the fringing over both edges
    bz, nu, nv = (0.022, 0.078, 0.0, 0.008), 90, 40
    acc = grid_idw(cen, mag, 0, 0.038, 0.052, 1, 2, nu, nv, bz, 0.0035, -25.0)
    write_grid("field_xz.dat", acc, bz, nu, nv, 1000.0, -25.0)
    nrm = [[(m + 25.0) / 25.0 for m in r] for r in acc]
    write_contours("cont_xz.dat", contours(nrm, bz, nu, nv, 1000.0, lv))


if __name__ == "__main__":
    if not os.path.exists(BIN):
        sys.exit("build the release binary first")
    conditioning()
    sweep()
    fields()
    print(mesh_cut(), 'mesh polygons')
    print("wrote cond.dat, patch.dat, field and contour grids")
