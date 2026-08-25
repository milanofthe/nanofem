#!/usr/bin/env python3
"""Regenerate the report data files from the release binary.

Run from the repository root after `cargo build --release`:

    python3 report/data/gen.py

Writes cond.dat, patch.dat, mesh3d.tikz and slice3d.tikz next to this file.
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
SUB = ("mat sub eps 2.2 tand 0.001\nmat sub_pml eps 2.2 tand 0.001\n"
       "pml sub_pml 3 3 3\npml air_pml 3 3 3\npec pec\n"
       "port 1 feed 0 0 1 200\n")

# The 3d figures are one scene seen from one place: same direction, same
# millimetres per millimetre, same origin, so the mesh and the field cut
# overlay. Geometry in mm, as in models/patch.geo.
H, SUBH, PMLT = 1.57, 0.00157, 10.0
DOM = (0.0, 90.0, 0.0, 100.0, 0.0, 31.57)
AZ, EL, WIDTH = 35.0, 20.0, 138.0


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
    plain = ["1.3e9", "1.5e8", "1.3e7", "1.5e6", "1.3e5", "1.5e4",
             "1.3e3", "3.0e2", "3.8e3", "2.0e5"]
    with open(os.path.join(HERE, "cond.dat"), "w") as f:
        f.write("f cond plain\n")
        for (fr, c), p in zip(rows, plain):
            f.write("%s %s %s\n" % (fr, c, p))


def sweep():
    """Reflection of the reference antenna, returning where it resonates."""
    out = run("mesh %s\n%ssweep lin 2.1e9 2.8e9 36\n" % (MESH, SUB))
    rows = []
    for line in out.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        v = line.split()
        m = math.hypot(float(v[1]), float(v[2]))
        rows.append((float(v[0]), 20 * math.log10(max(m, 1e-12))))
    with open(os.path.join(HERE, "patch.dat"), "w") as f:
        f.write("f db\n")
        for fr, db in rows:
            f.write("%.6e %.4f\n" % (fr, db))
    return min(rows, key=lambda r: r[1])[0]


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


def read_msh():
    """Nodes and tagged tetrahedra of the reference mesh, in mm."""
    ls = open(MESH).read().splitlines()
    i = ls.index("$Nodes")
    pos = {}
    for k in range(int(ls[i + 1])):
        t = ls[i + 2 + k].split()
        pos[int(t[0])] = tuple(float(x) * 1000.0 for x in t[1:4])
    i = ls.index("$Elements")
    tets = []
    for k in range(int(ls[i + 1])):
        t = ls[i + 2 + k].split()
        if t[1] == "4":
            nt = int(t[2])
            tets.append(([int(x) for x in t[3 + nt:]], int(t[3])))
    return pos, tets


def view():
    """The shared projection, as a depth key and a point formatter.

    The scale comes from the full domain, not from what a figure happens to
    draw, so every figure that uses it is at the same size and the same
    origin as the others.
    """
    a, e = math.radians(AZ), math.radians(EL)
    ex = (-math.sin(a), math.cos(a), 0.0)
    ey = (-math.cos(a) * math.sin(e), -math.sin(a) * math.sin(e), math.cos(e))
    ez = (math.cos(a) * math.cos(e), math.sin(a) * math.cos(e), math.sin(e))
    dot = lambda p, q: p[0] * q[0] + p[1] * q[1] + p[2] * q[2]
    x0, x1, y0, y1, z0, z1 = DOM
    pj = [(dot(c, ex), dot(c, ey)) for c in
          [(x, y, z) for x in (x0, x1) for y in (y0, y1) for z in (z0, z1)]]
    sc = WIDTH / (max(q[0] for q in pj) - min(q[0] for q in pj))
    ox, oy = min(q[0] for q in pj), min(q[1] for q in pj)
    P = lambda p: ((dot(p, ex) - ox) * sc, (dot(p, ey) - oy) * sc)
    D = lambda p: "(%.2fmm,%.2fmm)" % P(p)
    return (lambda p: dot(p, ez)), D


def boxwire(D, box=DOM):
    """The twelve edges of a box, thin and light."""
    c = [(x, y, z) for x in box[0:2] for y in box[2:4] for z in box[4:6]]
    return ["\\draw[line width=0.25pt, black!40] %s -- %s;" % (D(c[i]), D(c[j]))
            for i in range(8) for j in range(i + 1, 8)
            if sum(1 for k in range(3) if c[i][k] != c[j][k]) == 1]


def write(name, lines):
    with open(os.path.join(HERE, name), "w") as f:
        f.write("% generated by report/data/gen.py, do not edit\n"
                + "\n".join(lines) + "\n")
    return len(lines)


def mesh3d(name="mesh3d.tikz", cut=50.0):
    """Cutaway render of the mesh as painter sorted, filled polygons.

    Half the domain is removed, the boundary faces of what remains are
    projected, sorted back to front and written as tikz paths. Filling each
    face opaque is what makes the interior read as solid rather than as a
    wireframe. The domain is outlined first, so what the cut took away stays
    visible as an empty box.
    """
    pos, tets = read_msh()
    depth, D = view()
    keep = [(t, g) for t, g in tets if sum(pos[v][1] for v in t) / 4.0 <= cut]
    faces = {}
    for t, g in keep:
        for f in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
            k = tuple(sorted(t[i] for i in f))
            if k in faces:
                faces[k] = None
            else:
                faces[k] = g
    proj = []
    for k, g in faces.items():
        if g is None:
            continue
        p3 = [pos[v] for v in k]
        proj.append((sum(depth(q) for q in p3) / 3.0, p3, g))
    proj.sort(key=lambda r: r[0])
    L = boxwire(D)
    for _, p3, g in proj:
        # the two absorbing groups are shaded, the physical region is white
        fill = "black!12" if g in (3, 4) else "white"
        L.append("\\filldraw[fill=%s, draw=black, line width=0.15pt] %s -- cycle;"
                 % (fill, " -- ".join(D(q) for q in p3)))
    write(name, L)
    return len(proj)


def slice3d(segs, zcut, plane, patch, name="slice3d.tikz"):
    """One cut plane in the same scene, with its contours drawn on it.

    The contours are computed in the two dimensions of the plane itself, so
    the lines are exact; only the projection is shared. The plane is opaque
    and spans the interior, which is why the domain outline stands off it by
    the thickness of the absorbing layer.
    """
    _, D = view()
    x0, x1, y0, y1 = plane
    pxa, pxb, pya, pyb, ph = patch
    L = boxwire(D)
    L.append("\\filldraw[fill=white, draw=black, line width=0.3pt] %s -- cycle;"
             % " -- ".join(D(q) for q in [(x0, y0, zcut), (x1, y0, zcut),
                                          (x1, y1, zcut), (x0, y1, zcut)]))
    for p, q in segs:
        L.append("\\draw[line width=0.3pt] %s -- %s;"
                 % (D((p[0], p[1], zcut)), D((q[0], q[1], zcut))))
    L.append("\\draw[line width=0.8pt] %s -- cycle;"
             % " -- ".join(D(q) for q in [(pxa, pya, ph), (pxb, pya, ph),
                                          (pxb, pyb, ph), (pxa, pyb, ph)]))
    write(name, L)
    return len(segs)


def fields(freq):
    """A cut of |E| through the middle of the substrate at the resonance."""
    vtk = os.path.join(HERE, "_field.vtk")
    run("mesh %s\n%ssweep lin %.6e %.6e 1\nfield %s %.6e\n" % (MESH, SUB, freq, freq, vtk, freq))
    cen, mag = read_vtk(vtk)
    os.remove(vtk)
    # the interior of the domain, that is everything the absorbing ring leaves
    t, sx, sy = PMLT / 1000.0, 0.090, 0.100
    plane = (PMLT, 90.0 - PMLT, PMLT, 100.0 - PMLT)
    bx, nu, nv = (t, sx - t, t, sy - t), 90, 90
    acc = grid_idw(cen, mag, 2, 0.0, SUBH, 0, 1, nu, nv, bx, 0.005, -20.0)
    segs = contours([[(m + 20.0) / 20.0 for m in r] for r in acc], bx, nu, nv,
                    1000.0, [0.35, 0.55, 0.75, 0.9])
    return slice3d(segs, H / 2.0, plane, (20.8, 69.2, 31.6, 68.4, H))


if __name__ == "__main__":
    if not os.path.exists(BIN):
        sys.exit("build first: cargo build --release")
    conditioning()
    fr = sweep()
    print("%d 3d faces, %d contour segments, resonance %.3f GHz"
          % (mesh3d(), fields(fr), fr / 1e9))
    print("wrote cond.dat, patch.dat, mesh3d.tikz and slice3d.tikz")
