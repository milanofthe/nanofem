# nanofem

A 3D finite element electromagnetic field solver in one Rust source file,
capped at 1000 lines of code. No dependencies, std only. A test counts the
nonblank, noncomment lines of src/main.rs and fails above 1000. Current
count: 954.

The repository is educational. The report in report/ derives every algorithm
in the solver, explains the design decisions the budget forced, and maps
both to the code section by section.

## Build and run

    cargo build --release
    gmsh -3 -format msh22 -o models/patch.msh models/patch.geo
    target/release/nanofem models/patch.nfm

Output is Touchstone on stdout. Two models are in models/, each a Gmsh .geo
with a matching deck: an edge fed 2.45 GHz patch antenna terminated by a
PML, 37676 tetrahedra at 3.7 s per frequency, and a shielded line with a
lumped port at each end, 22775 tetrahedra at 8.1 s for 21 frequencies on
eight threads. Both exercise the solver, neither is a converged design.

## Scope

nanofem solves the time harmonic curl-curl equation for the electric field
with first order Nedelec edge elements on tetrahedra and reports scattering
parameters at lumped ports. One card per line, `*` starts a comment. Group
names refer to the physical groups defined in the mesh.

| Card | Meaning |
|---|---|
| `mesh <path>` | Gmsh .msh v2.2 ASCII mesh, path relative to the deck |
| `mat <group> eps <er> [tand <d>] [mur <mr>] [sigma <s>]` | material of a volume group, sigma in S/m |
| `pec <group> ...` | perfect electric conductor surfaces |
| `abc <group> ...` | first order absorbing boundary surfaces |
| `metal <group> <sigma>` | lossy conductor sheet, conductivity in S/m |
| `pml <group> <ax> <ay> <az>` | PML volume, imaginary coordinate stretch per axis |
| `port <n> <group> <jx> <jy> <jz> <z0>` | lumped port: number, surface group, voltage direction, reference impedance |
| `sweep lin <f0> <f1> <npoints>` | frequency sweep in Hz |
| `output <s\|z\|y\|lq>` | what to print, default s |
| `field <path.vtk> <f>` | E field snapshot at f with port 1 driven, legacy VTK |

Volume groups without a `mat` card are vacuum. Surfaces without a role are
natural boundaries, which for the curl-curl equation means PMC, so a
magnetic symmetry plane needs no card. Ports are rectangular sheets: the
direction vector points from one terminal to the other, the height is the
mesh extent along it, the width follows from the face area. Put nothing or
PEC behind a PML. `z` and `y` print the impedance and admittance matrices,
`lq` reads each port as a coil, all three as comma separated values.
Triangles, tetrahedra and their physical groups are read, the rest of the
mesh file is ignored. Not supported: modal waveguide ports, adaptive
refinement, elements beyond first order, dispersive materials.

## Algorithms

Six Whitney edge functions per tetrahedron. Element vertices are sorted on
input, so every local edge runs from the lower to the higher global node and
no orientation sign appears anywhere in the solver. PEC surfaces eliminate
their unknowns, the absorbing boundary and the port sheets enter as face
mass matrices, and a PML is a complex coordinate stretch applied as a
diagonal tensor to eps and mu. Loss stays frequency independent in the
assembly: a loss tangent scales the imaginary part of eps, a volume
conductivity adds j k0 eta0 sigma, and a conductor sheet collapses to a
constant times the square root of k0, so every matrix entry carries four
coefficients against the basis 1, k0, k0 squared and the square root of k0
and the assembly runs once for an entire sweep. The system is complex
symmetric and solved directly: geometric nested dissection, then a sparse
LDL^T, equilibrated with the inverse square root of the diagonal, which cuts
the pivot spread by more than a factor of six at low frequency, and one step
of iterative refinement per solve against the unscaled matrix, which
produces a measured residual. Frequencies run in parallel, the thread count
capped by a memory budget since each thread holds its own factorization.
Derivations are in report/nanofem.pdf.

Diagnostics go to stderr: how every physical group was understood, including
the ones the deck never names, then the worst pivot spread, a free lower
bound on the condition number, and the worst relative residual. A malformed
mesh is reported rather than indexed past the end; the deck is otherwise not
checked.

## Tests

`cargo test` runs 12 integration tests in tests/cli.rs against the built
binary, almost all against closed form results: a matched parallel plate TEM
line in magnitude and phase, a lossy dielectric line, a shorted line, a
mismatched port, an absorbing wall and a PML slab each terminating a TEM
wave, a PEC cavity at its analytic mode, a conductive filling, lossy plates,
a shorted line read as a coil, a malformed mesh producing a diagnostic
rather than a backtrace, and the LOC budget guard. Tests and comments do not
count toward the budget.

## Report

    tectonic report/nanofem.tex

The prebuilt PDF is committed at report/nanofem.pdf. The figures read data
files from report/data/, regenerated from the release binary with
python3 report/data/gen.py.

## License

MIT.
