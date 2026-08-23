# nanofem

A headless 3D finite element electromagnetic field solver in one Rust source
file, capped at 1000 lines of code. No dependencies, std only. A test counts
the nonblank, noncomment lines of src/main.rs and fails above 1000. Current
count: 954.

The repository is educational. The report in report/ derives the
formulation, explains every design decision the budget forced, and maps both
to the code section by section.

nanofem solves the time harmonic curl-curl equation for the electric field
with first order Nedelec edge elements on tetrahedra and reports scattering
parameters at lumped ports. The mesh comes from Gmsh, the setup is a small
text deck that maps physical group names to materials, boundaries and ports.

## Build and run

    cargo build --release
    target/release/nanofem antenna.nfm

Output goes to stdout, diagnostics to stderr. The default output is
Touchstone. The `output` card switches to derived port quantities as comma
separated values with a header: `z` and `y` are the impedance and admittance
matrices, `lq` reads each port as a coil and gives its inductance and
quality factor.

## Deck format

One card per line, `*` starts a comment. Group names refer to the physical
groups defined in the Gmsh mesh.

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
direction vector points from one terminal to the other, the port height is
the mesh extent along that direction, and the width follows from the face
area. Put nothing or PEC behind a PML.

The mesh must be Gmsh version 2.2 ASCII (`gmsh -3 -format msh22 model.geo`).
Triangles, tetrahedra and their physical groups are read, everything else in
the file is ignored.

Not supported: modal waveguide ports, adaptive refinement, elements beyond
first order, dispersive materials.

## Algorithms

The electric field is discretized with the six Whitney edge functions per
tetrahedron. PEC surfaces eliminate their unknowns, a first order absorbing
boundary and the port sheets enter as face mass matrices, and a PML region
is a complex coordinate stretch applied as a diagonal tensor to both eps and
mu, which keeps the wave impedance matched while the field decays.

Element vertices are sorted on input. Every local edge then runs from the
lower to the higher global node, so neighboring elements agree on their
Whitney functions by construction and no orientation signs appear anywhere
in the solver.

Loss comes in three forms and all of them stay frequency independent in the
assembly. A dielectric loss tangent scales the imaginary part of eps. A
volume conductivity adds a conduction current, which works out to exactly
+j k0 eta0 sigma. A `metal` surface carries the impedance of a good
conductor, whose surface resistance grows with the square root of frequency,
and j k0 eta0 / Zs collapses to a constant times the square root of k0. Each
matrix entry therefore keeps four coefficients against the basis 1, k0, k0
squared and the square root of k0, and the whole assembly happens once for
an entire sweep.

The system is complex symmetric and is solved directly: a geometric nested
dissection ordering, then a sparse LDL^T. Before factoring, the matrix is
equilibrated symmetrically with the inverse square root of its diagonal,
which cuts the pivot spread by more than a factor of six at low frequency.
Every solve is followed by one step of iterative refinement against the
unscaled matrix, which both repairs the accuracy a weak pivot costs and
produces a measured residual. The frequencies of a sweep are independent and
run in parallel, with the thread count capped by a memory budget because
each thread holds its own factorization.

Derivations are in report/nanofem.pdf.

## Diagnostics

Before solving, nanofem prints how it understood every physical group,
including the ones the deck never names. A group the deck does not mention
becomes vacuum or a natural PMC wall, which is valid input and therefore
not reachable by any check.

After the sweep it prints the worst pivot spread and the worst relative
residual. The pivot spread is a free lower bound on the condition number: it
grows like one over frequency squared towards low frequency, where the
curl-curl operator loses the mass term that regularizes its nullspace, and
it grows again once the mesh gets coarse against the wavelength. The
residual states how well the system was solved, measured rather than
assumed.

A malformed mesh file is reported rather than indexed past the end. The
deck is otherwise not checked: whether a model makes sense is decided by
whatever generates it.

## Models

models/ holds an edge fed 2.45 GHz microstrip patch antenna. patch.geo
builds the geometry and mesh with Gmsh, patch.nfm is the matching deck.

    gmsh -3 -format msh22 -o models/patch.msh models/patch.geo
    target/release/nanofem models/patch.nfm

The sweep shows the resonance as an S11 dip at 2.45 GHz on the default mesh.
That number is not converged: refining the mesh moves it to 2.50 GHz, and
part of the remaining drift is the absorbing boundary, which sits only about
a third of a wavelength above the patch. The model exercises the solver, it
is not a converged antenna design.

## Tests

`cargo test` runs 12 integration tests in tests/cli.rs against the built
binary, almost all against closed form results: a matched parallel plate TEM
line in magnitude and phase, a lossy dielectric filled line against the
analytic attenuation, a shorted line reflecting with unit magnitude and the
right phase, a deliberately mismatched port against the impedance
transformation, an absorbing wall and a PML slab each terminating a TEM
wave, a PEC box cavity at its analytic mode frequency, a conductive filling
against the analytic propagation constant, lossy plates against the analytic
surface resistance, a shorted line whose impedance, admittance and
inductance match the closed form for a coil, a malformed mesh producing a
diagnostic rather than a backtrace, and the LOC budget guard. Tests and
comments do not count toward the budget.

## Report

    tectonic report/nanofem.tex

The prebuilt PDF is committed at report/nanofem.pdf.

## License

MIT.
