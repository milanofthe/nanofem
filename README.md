# nanofem

A headless 3D finite element electromagnetic field solver in one Rust source
file, capped at 1000 lines of code. No dependencies, std only. A test counts
the nonblank, noncomment lines of src/main.rs and fails above 1000. The
budget covers the solver alone: tests, comments and the models in models/ do
not count toward it.

nanofem solves the time harmonic curl-curl equation for the electric field
with first order Nedelec (Whitney) edge elements on tetrahedra and computes
scattering parameters at lumped ports. The mesh comes from Gmsh, the solver
setup is a small text deck that maps physical group names to materials,
boundaries and ports.

## Build and run

    cargo build --release
    target/release/nanofem antenna.nfm

Output is Touchstone data on stdout, diagnostics go to stderr. The `output`
card switches to derived port quantities, printed as comma separated values
with a header: `z` and `y` give the impedance and admittance matrices, `lq`
reads each port as a coil and gives its inductance and quality factor. Those
are the numbers an inductor or interconnect extraction actually wants.

## Deck format

One card per line, `*` starts a comment. Physical group names refer to the
names defined in the Gmsh mesh.

| Card | Meaning |
|---|---|
| `mesh <path>` | Gmsh .msh v2.2 ASCII mesh, path relative to the deck |
| `mat <group> eps <er> [tand <d>] [mur <mr>] [sigma <s>]` | material of a volume group, sigma in S/m |
| `pec <group> ...` | perfect electric conductor surfaces |
| `abc <group> ...` | first order absorbing boundary surfaces |
| `pml <group> <ax> <ay> <az>` | PML volume: imaginary coordinate stretch per axis |
| `metal <group> <sigma>` | lossy conductor sheet of conductivity sigma in S/m |
| `port <n> <group> <jx> <jy> <jz> <z0>` | lumped port: number, surface group, voltage direction, reference impedance |
| `sweep lin <f0> <f1> <npoints>` | frequency sweep in Hz |
| `field <path.vtk> <f>` | E field snapshot at f with port 1 driven, legacy VTK |
| `output <s\|z\|y\|lq>` | what to print, default s |

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
mesh. That number is not converged: refining the mesh moves it to 2.48 GHz,
and part of the remaining drift is the first order absorbing boundary,
which sits only about a third of a wavelength above the patch. The model is
meant to exercise the solver, not to be a converged antenna design.

## Mesh

Gmsh version 2.2 ASCII format (`gmsh -3 -format msh22 model.geo`). The
solver reads triangles and tetrahedra and their physical groups; everything
else in the file is ignored.

## Scope

PEC, first order ABC, PML regions, natural PMC, lossy dielectrics and
magnetics per region, conductive volumes, lossy conductor sheets with the
skin effect, lumped rectangular ports with S, Z, Y and LQ extraction, E
field export to VTK for ParaView.

Loss comes in three forms. A dielectric loss tangent scales the imaginary
part of eps. A volume conductivity adds a conduction current, which enters
the system as a term linear in the wave number. A `metal` surface carries
the impedance of a good conductor, whose surface resistance grows with the
square root of frequency; every matrix entry therefore keeps four
frequency independent coefficients, against 1, k0, k0 squared and the
square root of k0, so the whole assembly still happens once for an entire
sweep.

A PML region is declared in the deck by naming a mesh volume and the axes
along which it absorbs; put nothing or PEC behind it. The system is solved
directly with a complex symmetric sparse LDL^T after a geometric nested
dissection ordering, and the frequencies of a sweep run in parallel across
threads.

No modal waveguide ports, no adaptive refinement, first order elements
only.

## Diagnostics

Before solving, nanofem prints on stderr how it understood every physical
group in the mesh, including the ones the deck never names, since a group
that silently defaults to vacuum or to a natural PMC wall is the one
mistake input validation cannot catch. After the sweep it prints the worst
pivot spread of the factorizations, a free lower bound on the condition
number. That number grows like one over frequency squared towards low
frequency, where the curl curl operator loses the mass term that
regularizes its nullspace, and it grows again once the mesh gets coarse
against the wavelength, so it tells you whether the run was in the
trustworthy regime.

A malformed mesh file is reported rather than indexed past the end, but
nanofem does not otherwise police the deck. Checking that a model makes
sense is the job of whatever generates the deck, the same division of
labour nanospice uses for netlists.

## Validation

The integration tests build structured meshes and check against closed form
results: a matched parallel plate TEM line (S11, S21 magnitude and phase), a
lossy dielectric filled line against the analytic attenuation and phase, a
shorted line reflecting with unit magnitude and the right phase, a
deliberately mismatched port against the impedance transformation, an
absorbing wall terminating a TEM wave, a PML slab doing the same, a PEC box
cavity resonating at the analytic mode frequency, a conductive filling
attenuating by the analytic propagation constant, lossy plates attenuating
by the analytic surface resistance, and a shorted line whose impedance,
admittance and inductance match the closed form for a coil.
