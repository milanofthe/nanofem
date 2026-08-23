// Edge fed microstrip patch antenna, 2.45 GHz design on eps_r 2.2, h 1.57 mm.
// Grounded substrate with an air box above, both wrapped in a perfectly
// matched layer: a ring of thickness t around the footprint and a slab of the
// same thickness on top. The feed is a vertical lumped port sheet from the
// ground plane to the patch edge. A PML needs no standoff from the radiator,
// so the air region is a third of the height an absorbing boundary would
// need. All dimensions in mm, scaled to meters on output.
// Mesh with: gmsh -3 -format msh22 -o patch.msh patch.geo

SetFactory("OpenCASCADE");

sx = 90; sy = 100;     // footprint
h  = 1.57;             // substrate height
Ha = 14;               // air above the substrate
t  = 8;                // pml thickness
pw = 48.4; pl = 36.8;  // patch width (x) and resonant length (y)
fw = 4;                // feed width
H  = h + Ha + t;       // total height

// substrate, split into the interior and the absorbing ring
Box(1) = {0, 0, 0, sx, sy, h};
Box(2) = {t, t, 0, sx - 2*t, sy - 2*t, h};
// air, split the same way; the outer volume also carries the top slab
Box(3) = {0, 0, h, sx, sy, H - h};
Box(4) = {t, t, h, sx - 2*t, sy - 2*t, Ha};

px = (sx - pw) / 2; py = (sy - pl) / 2;
Rectangle(100) = {px, py, h, pw, pl};
fx = sx / 2 - fw / 2;
Rectangle(101) = {fx, py, 0, fw, h};
Rotate{ {1, 0, 0}, {fx, py, 0}, Pi/2 }{ Surface{101}; }
BooleanFragments{ Volume{1, 2, 3, 4}; Delete; }{ Surface{100, 101}; Delete; }

eps = 1e-3;
inner()  = Volume In BoundingBox{t-eps, t-eps, -eps,
                                 sx-t+eps, sy-t+eps, h+Ha+eps};
subv()   = Volume In BoundingBox{-eps, -eps, -eps, sx+eps, sy+eps, h+eps};
airv()   = Volume In BoundingBox{-eps, -eps, h-eps, sx+eps, sy+eps, H+eps};
subin()  = Volume In BoundingBox{t-eps, t-eps, -eps, sx-t+eps, sy-t+eps, h+eps};
airin()  = Volume In BoundingBox{t-eps, t-eps, h-eps,
                                 sx-t+eps, sy-t+eps, h+Ha+eps};

patch()  = Surface In BoundingBox{px-eps, py-eps, h-eps,
                                  px+pw+eps, py+pl+eps, h+eps};
feed()   = Surface In BoundingBox{fx-eps, py-eps, -eps, fx+fw+eps, py+eps, h+eps};
gnd()    = Surface In BoundingBox{-eps, -eps, -eps, sx+eps, sy+eps, eps};

// The absorbing regions are what remains of each layer once the interior is
// taken out. Gmsh has no set difference on volume lists, so subtract by hand.
subring() = {};
For i In {0 : #subv()-1}
  keep = 1;
  For j In {0 : #subin()-1}
    If (subv(i) == subin(j))
      keep = 0;
    EndIf
  EndFor
  If (keep == 1)
    subring() += subv(i);
  EndIf
EndFor
airring() = {};
For i In {0 : #airv()-1}
  keep = 1;
  For j In {0 : #airin()-1}
    If (airv(i) == airin(j))
      keep = 0;
    EndIf
  EndFor
  If (keep == 1)
    airring() += airv(i);
  EndIf
EndFor

Physical Volume("sub") = {subin()};
Physical Volume("air") = {airin()};
Physical Volume("sub_pml") = {subring()};
Physical Volume("air_pml") = {airring()};
Physical Surface("pec") = {patch(), gnd()};
Physical Surface("feed") = {feed()};

// Grade the element size with distance from the metal: fine on the patch and
// the feed, growing smoothly outwards. A uniform maximum size instead jumps
// from the substrate thickness straight to the coarse size and produces
// slivers spanning the air region.
Field[1] = Distance;
Field[1].SurfacesList = {patch(), feed()};
Field[1].Sampling = 80;
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 3.0;
Field[2].SizeMax = 9.0;
Field[2].DistMin = 3.0;
Field[2].DistMax = 26.0;
Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.Optimize = 1;
Mesh.ScalingFactor = 0.001;
