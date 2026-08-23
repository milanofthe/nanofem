#!/usr/bin/env python3
"""Regenerate the report data files from the release binary.

Run from the repository root after `cargo build --release`:

    python3 report/data/gen.py

Writes cond.dat, patch.dat, field_xy.dat and field_xz.dat next to this file.
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
    plain = ["1.0e9", "1.2e8", "1.0e7", "1.2e6", "1.0e5", "1.2e4",
             "1.0e3", "5.0e2", "3.2e4", "1.8e5"]
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


def grid_nn(cen, mag, axis, lo, hi, u, v, nu, nv, bounds, floor_db=-30.0):
    """Nearest neighbour sample of |E| in dB onto a regular nu by nv grid.

    The mesh is coarser than any useful raster, so averaging leaves holes.
    Nearest neighbour gives a piecewise constant map with no gaps, which also
    suits the stylized look of the figures. The scale is dB relative to the
    peak of the cut, clipped at floor_db, because the fringing field above a
    patch is two orders of magnitude below the field inside the substrate.
    """
    u0, u1, v0, v1 = bounds
    pts = [(c[u], c[v], m) for c, m in zip(cen, mag) if lo <= c[axis] <= hi]
    if not pts:
        sys.exit("empty slab")
    peak = max(m for _, _, m in pts) or 1.0
    out = []
    for r in range(nv):
        y = v0 + (v1 - v0) * (r + 0.5) / nv
        row = []
        for k in range(nu):
            x = u0 + (u1 - u0) * (k + 0.5) / nu
            best, bd = 0.0, float("inf")
            for cu, cv, m in pts:
                d = (cu - x) ** 2 + (cv - y) ** 2
                if d < bd:
                    bd, best = d, m
            db = 20.0 * math.log10(max(best / peak, 10 ** (floor_db / 20.0)))
            row.append(max(db, floor_db))
        out.append(row)
    return out


def write_grid(name, acc, bounds, nu, nv, scale, floor_db=-30.0):
    u0, u1, v0, v1 = bounds
    with open(os.path.join(HERE, name), "w") as f:
        f.write("x y e\n")
        for r in range(nv):
            for k in range(nu):
                x = (u0 + (u1 - u0) * (k + 0.5) / nu) * scale
                y = (v0 + (v1 - v0) * (r + 0.5) / nv) * scale
                f.write("%.4f %.4f %.4f\n" % (x, y, (acc[r][k] - floor_db) / (-floor_db)))
            f.write("\n")


def fields():
    """Two cuts of |E| through the antenna at its resonance."""
    vtk = os.path.join(HERE, "_field.vtk")
    run("mesh %s\n%ssweep lin 2.4e9 2.4e9 1\nfield %s 2.4e9\n" % (MESH, SUB, vtk))
    cen, mag = read_vtk(vtk)
    os.remove(vtk)
    h = 0.00157
    # inside the substrate, looking down on the patch, cropped to the region
    # that carries field
    bx = (0.012, 0.078, 0.022, 0.078)
    acc = grid_nn(cen, mag, 2, 0.0, h, 0, 1, 44, 38, bx, -20.0)
    write_grid("field_xy.dat", acc, bx, 44, 38, 1000.0, -20.0)
    # vertical cut along the resonant direction, through the middle of the
    # patch: shows the half wave in the substrate and the fringing over both
    # radiating edges
    bz = (0.022, 0.078, 0.0, 0.006)
    acc = grid_nn(cen, mag, 0, 0.040, 0.050, 1, 2, 44, 16, bz, -25.0)
    write_grid("field_xz.dat", acc, bz, 44, 16, 1000.0, -25.0)


if __name__ == "__main__":
    if not os.path.exists(BIN):
        sys.exit("build the release binary first")
    conditioning()
    sweep()
    fields()
    print("wrote cond.dat, patch.dat, field_xy.dat, field_xz.dat")
