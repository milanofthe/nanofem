// Shielded microstrip line on eps_r 2.2, h 1.57 mm. A strip over a ground
// plane inside a metal box, with a lumped port at each end reaching from the
// ground plane to the strip. W = 4.8 mm is the Hammerstad width for 50 ohm on
// this cross section, where the closed form gives eps_eff = 1.870 and
// Z0 = 50.5 ohm. That closed form is for an open microstrip and takes the
// port reference plane at the end face; neither holds here, so the sweep is
// not a check against it. All dimensions in mm, scaled to meters on output.
// Mesh with: gmsh -3 -format msh22 -o microstrip.msh microstrip.geo

SetFactory("OpenCASCADE");

lx = 44;               // line length
ly = 20;               // box width
h  = 1.57;             // substrate height
ha = 8;                // air above the substrate
w  = 4.8;              // strip width
H  = h + ha;

Box(1) = {0, 0, 0, lx, ly, h};
Box(2) = {0, 0, h, lx, ly, ha};
y0 = (ly - w) / 2;
Rectangle(100) = {0, y0, h, lx, w};
// port sheets: built flat, then swung into the y-z plane so that they span
// the strip width in y and the substrate height in z
Rectangle(101) = {0, y0, 0, h, w};
Rotate{ {0, 1, 0}, {0, y0, 0}, -Pi/2 }{ Surface{101}; }
Rectangle(102) = {0, y0, 0, h, w};
Rotate{ {0, 1, 0}, {0, y0, 0}, -Pi/2 }{ Surface{102}; }
Translate{ lx, 0, 0 }{ Surface{102}; }
BooleanFragments{ Volume{1, 2}; Delete; }{ Surface{100, 101, 102}; Delete; }

eps = 1e-3;
subv()  = Volume In BoundingBox{-eps, -eps, -eps, lx+eps, ly+eps, h+eps};
airv()  = Volume In BoundingBox{-eps, -eps, h-eps, lx+eps, ly+eps, H+eps};
strip() = Surface In BoundingBox{-eps, y0-eps, h-eps, lx+eps, y0+w+eps, h+eps};
p1()    = Surface In BoundingBox{-eps, y0-eps, -eps, eps, y0+w+eps, h+eps};
p2()    = Surface In BoundingBox{lx-eps, y0-eps, -eps, lx+eps, y0+w+eps, h+eps};
gnd()   = Surface In BoundingBox{-eps, -eps, -eps, lx+eps, ly+eps, eps};
top()   = Surface In BoundingBox{-eps, -eps, H-eps, lx+eps, ly+eps, H+eps};
ym()    = Surface In BoundingBox{-eps, -eps, -eps, lx+eps, eps, H+eps};
yp()    = Surface In BoundingBox{-eps, ly-eps, -eps, lx+eps, ly+eps, H+eps};

Physical Volume("sub") = {subv()};
Physical Volume("air") = {airv()};
Physical Surface("pec") = {strip(), gnd(), top(), ym(), yp()};
Physical Surface("p1") = {p1()};
Physical Surface("p2") = {p2()};

// Fine on the strip, coarser towards the shield, about three elements through
// the substrate. Refining from 1.4 mm to 0.45 mm moves the extracted phase
// constant by two percent, refining from 0.6 mm to 0.45 mm by 0.2 percent.
Field[1] = Distance;
Field[1].SurfacesList = {strip()};
Field[1].Sampling = 60;
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 0.6;
Field[2].SizeMax = 4.0;
Field[2].DistMin = 0.7;
Field[2].DistMax = 9.0;
Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.Optimize = 1;
Mesh.ScalingFactor = 0.001;
