# nanofem

A headless 3D finite element electromagnetic field solver in one Rust source
file, capped at 1000 lines of code. No dependencies, std only. A test counts
the nonblank, noncomment lines of src/main.rs and fails above 1000. The
budget covers the solver alone: tests, comments and the models in models/ do
not count toward it.

nanofem solves the time harmonic curl-curl equation for the electric field
with Nedelec edge elements on tetrahedra, first or second order, and
computes scattering parameters at lumped ports. The mesh comes from Gmsh,
the solver setup is a small text deck that maps physical group names to
materials, boundaries and ports.

## Build and run

    cargo build --release
    target/release/nanofem antenna.nfm

Output is Touchstone data on stdout, diagnostics go to stderr.

## Deck format

One card per line, `*` starts a comment. Physical group names refer to the
names defined in the Gmsh mesh.

| Card | Meaning |
|---|---|
| `mesh <path>` | Gmsh .msh v2.2 ASCII mesh, path relative to the deck |
| `mat <group> eps <er> [tand <d>] [mur <mr>]` | material of a volume group |
| `pec <group> ...` | perfect electric conductor surfaces |
| `abc <group> ...` | first order absorbing boundary surfaces |
| `pml <group> <ax> <ay> <az>` | PML volume: imaginary coordinate stretch per axis |
| `port <n> <group> <jx> <jy> <jz> <z0>` | lumped port: number, surface group, voltage direction, reference impedance |
| `sweep lin <f0> <f1> <npoints>` | frequency sweep in Hz |
| `order <1\|2>` | element order, default 1 |
| `field <path.vtk> <f>` | E field snapshot at f with port 1 driven, legacy VTK |

Volume groups without a `mat` card are vacuum. Surfaces without a role are
natural boundaries, which for the curl-curl equation means PMC. Ports are
rectangular sheets: the direction vector points from one terminal to the
other, the port height is the mesh extent along that direction and the width
follows from the face area.

## Models

models/ contains an edge fed 2.45 GHz microstrip patch antenna: patch.geo
builds the geometry and mesh with Gmsh, patch.nfm is the matching deck.

    gmsh -3 -format msh22 -o models/patch.msh models/patch.geo
    target/release/nanofem models/patch.nfm

The sweep shows the resonance as an S11 dip at 2.40 GHz on the default
mesh. That number is not converged: refining the mesh moves it to 2.48 GHz
and adding `order 2` to the deck moves it to 2.53 GHz on the coarse mesh
and 2.56 GHz on a medium one. Part of the remaining drift is the first
order absorbing boundary, which sits only about a third of a wavelength
above the patch. The model is meant to exercise the solver, not to be a
converged antenna design.

## Mesh

Gmsh version 2.2 ASCII format (`gmsh -3 -format msh22 model.geo`). The
solver reads triangles and tetrahedra and their physical groups; everything
else in the file is ignored.

## Scope

PEC, first order ABC, PML regions, natural PMC, lossy dielectrics and
magnetics per region, lumped rectangular ports with S-parameter extraction,
E field export to VTK for ParaView. Direct solution with a complex symmetric
sparse LDL^T after a geometric nested dissection ordering, frequencies run
in parallel across threads. A PML region is declared in the deck by naming
a mesh volume and the axes along which it absorbs; put nothing or PEC
behind it.

Elements are first or second order. Second order adds the curl free edge
gradients and two functions per face, 20 per tetrahedron instead of 6.
Element vertices are sorted on input, so local edges and faces inherit the
global numbering and no orientation signs are needed anywhere. Second order
costs roughly five times the unknowns of the same mesh, so watch memory:
each thread holds its own factorization.

No modal waveguide ports, no adaptive refinement, no higher order geometry.

## Validation

The integration tests build structured meshes and check against closed form
results: a matched parallel plate TEM line (S11, S21 magnitude and phase), a
lossy dielectric filled line against the analytic attenuation and phase, a
shorted line reflecting with unit magnitude and the right phase, a
deliberately mismatched port against the impedance transformation, an
absorbing wall terminating a TEM wave, a PML slab doing the same, and a PEC
box cavity resonating at the analytic mode frequency.

Element order is pinned by the convergence rate itself. The phase a wave
accumulates over the TEM line is known exactly, so the error of a run is
known exactly, and theory says it grows as k^(2p+1). Tripling the frequency
must therefore multiply the error by about 27 at first order and about 243
at second, which is what the tests require. That rate is a sharp check: a
wrong sign or a mismatched face function between neighboring tetrahedra
still gives plausible looking S-parameters, but it destroys the rate.
