// Edge fed microstrip patch antenna, 2.45 GHz design on eps_r 2.2, h 1.57 mm.
// Substrate slab with the ground plane as the domain floor, air box above,
// ABC on the remaining outer faces. The feed is a vertical lumped port sheet
// from ground to the patch edge. All dimensions in mm, scaled to meters on
// output. Mesh with: gmsh -3 -format msh22 -o patch.msh patch.geo

SetFactory("OpenCASCADE");

sx = 90; sy = 100;     // footprint
h  = 1.57;             // substrate height
H  = 35;               // domain height
pw = 48.4; pl = 36.8;  // patch width (x) and resonant length (y)
fw = 4;                // feed width

Box(1) = {0, 0, 0, sx, sy, h};
Box(2) = {0, 0, h, sx, sy, H - h};
px = (sx - pw) / 2; py = (sy - pl) / 2;
Rectangle(100) = {px, py, h, pw, pl};
fx = sx / 2 - fw / 2;
Rectangle(101) = {fx, py, 0, fw, h};
Rotate{ {1, 0, 0}, {fx, py, 0}, Pi/2 }{ Surface{101}; }
BooleanFragments{ Volume{1, 2}; Delete; }{ Surface{100, 101}; Delete; }

eps = 1e-3;
subv()  = Volume In BoundingBox{-eps, -eps, -eps, sx+eps, sy+eps, h+eps};
airv()  = Volume In BoundingBox{-eps, -eps, h-eps, sx+eps, sy+eps, H+eps};
patch() = Surface In BoundingBox{px-eps, py-eps, h-eps, px+pw+eps, py+pl+eps, h+eps};
feed()  = Surface In BoundingBox{fx-eps, py-eps, -eps, fx+fw+eps, py+eps, h+eps};
gnd()   = Surface In BoundingBox{-eps, -eps, -eps, sx+eps, sy+eps, eps};
top()   = Surface In BoundingBox{-eps, -eps, H-eps, sx+eps, sy+eps, H+eps};
xm()    = Surface In BoundingBox{-eps, -eps, -eps, eps, sy+eps, H+eps};
xp()    = Surface In BoundingBox{sx-eps, -eps, -eps, sx+eps, sy+eps, H+eps};
ym()    = Surface In BoundingBox{-eps, -eps, -eps, sx+eps, eps, H+eps};
yp()    = Surface In BoundingBox{-eps, sy-eps, -eps, sx+eps, sy+eps, H+eps};

Physical Volume("sub") = {subv()};
Physical Volume("air") = {airv()};
Physical Surface("pec") = {patch(), gnd()};
Physical Surface("feed") = {feed()};
Physical Surface("open") = {top(), xm(), xp(), ym(), yp()};

// Grade the element size with distance from the metal: fine on the patch
// and the feed, growing smoothly into the air box. Without this the mesh
// jumps from the substrate thickness straight to the maximum size and
// produces slivers.
Field[1] = Distance;
Field[1].SurfacesList = {patch(), feed()};
Field[1].Sampling = 80;
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 3.0;
Field[2].SizeMax = 11.0;
Field[2].DistMin = 3.0;
Field[2].DistMax = 28.0;
Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.Optimize = 1;
Mesh.ScalingFactor = 0.001;
