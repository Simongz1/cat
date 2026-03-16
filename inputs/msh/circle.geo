// Mesh.RecombinationAlgorith = 1;
// Mesh.RecombineAll = 1;

SetFactory("OpenCASCADE");

radius = r;

w = 4 * radius;
h = w;

Rectangle(1) = {0, 0, 0, h, w, 0};
Rectangle(2) = {-2*w/30, -2*w/30, 0, -w/30, -w/30, 0};
Disk(3) = {w/2, h/2, 0, radius, radius};

BooleanDifference(4) = { Surface{1}; Delete; }{ Surface{3}; Delete; };

out1[] = Extrude {0, 0, radius/10} { Surface{4}; };
out2[] = Extrude {0, 0, radius/10} { Surface{2}; };

Physical Volume("vol") = { out1[1] };
Physical Volume("dead") = { out2[1] };

Field[1] = Box;
Field[1].Thickness = 100;
Field[1].VIn = r/15;
Field[1].VOut = 50;
Field[1].XMax = 100;
Field[1].YMax = 100;
Field[1].ZMax = 100;
Background Field = 1;

Physical Surface("front", 32) = {4, 10};
Physical Surface("right", 33) = {7};
Physical Surface("left", 34) = {6};
Physical Surface("top", 35) = {8};
Physical Surface("bottom", 36) = {5};
