# nanofem

A headless 3D finite element electromagnetic field solver in one Rust source
file, capped at 1000 lines of code. No dependencies, std only. A test counts
the nonblank, noncomment lines of src/main.rs and fails above 1000.

nanofem solves the time harmonic curl-curl equation for the electric field
with first order Nedelec (Whitney) edge elements on tetrahedra and computes
scattering parameters at lumped ports. The mesh comes from Gmsh, the solver
setup is a small text deck that maps physical group names to materials,
boundaries and ports.

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
| `port <n> <group> <jx> <jy> <jz> <z0>` | lumped port: number, surface group, voltage direction, reference impedance |
| `sweep lin <f0> <f1> <npoints>` | frequency sweep in Hz |

Volume groups without a `mat` card are vacuum. Surfaces without a role are
natural boundaries, which for the curl-curl equation means PMC. Ports are
rectangular sheets: the direction vector points from one terminal to the
other, the port height is the mesh extent along that direction and the width
follows from the face area.

## Example

examples/ contains an edge fed 2.45 GHz microstrip patch antenna: patch.geo
builds the geometry and mesh with Gmsh, patch.nfm is the matching deck.

    gmsh -3 -format msh22 -o examples/patch.msh examples/patch.geo
    target/release/nanofem examples/patch.nfm

The sweep shows the resonance as an S11 dip at 2.40 GHz on the default
mesh, about 2 percent below the mesh converged 2.45 GHz.

## Mesh

Gmsh version 2.2 ASCII format (`gmsh -3 -format msh22 model.geo`). The
solver reads triangles and tetrahedra and their physical groups; everything
else in the file is ignored.

## Scope

PEC, first order ABC, natural PMC, lossy dielectrics and magnetics per
region, lumped rectangular ports with S-parameter extraction. Direct
solution with a complex symmetric skyline LDL^T after reverse Cuthill-McKee
ordering, so the practical problem size is a few tens of thousands of
unknowns. No modal waveguide ports, no PML, no adaptive refinement, no field
export.

## Validation

The integration tests build structured meshes and check against closed form
results: a matched parallel plate TEM line (S11, S21 magnitude and phase), a
lossy dielectric filled line against the analytic attenuation and phase, a
shorted line reflecting with unit magnitude and the right phase, a
deliberately mismatched port against the impedance transformation, an
absorbing wall terminating a TEM wave, and a PEC box cavity resonating at
the analytic mode frequency.
