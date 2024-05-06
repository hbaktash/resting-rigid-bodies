// Deltoidal Icositetrahedron
// 4 sided faces = 24
// degree 4 vertices = 18
// degree 3 vertices = 8
// "halves" = 6

txt_depth = .09;
txt_size = .35;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 24; // scale factor which sets the diameter (distance from one [deg 4]-vertex to the vertex opposite) 
minko = 0.2; //chamfer the edges [0 = disabled]
roll = 0; // another effet to round the dice (intersects with a smaller sphere) [0 = disabled]
minkfn = 80; //sets the $fn variable for chamfer and the sphere

C0 = 0.773459080339013578400241246316;
C1 = 1.41421356237309504880168872421;

ordi = C1+minko; //original diameter
scafa = diam*.5/ordi; //scale factor


//labels = ["1","2","3","4","5","6.","7","8","9.","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24"]; // default ordering

// this solid can be but in halves in 3 different ways (yielding 6 halves)

//these labeling are compromises: halves are balanced; the opposed faces deviate by 3 points (for 6 pairs) and 1 point (for the other 6); 14 of the degree 4 vertices are balanced while the remaining 4 deviate by 1 point; 4 of the degree 3 vertices deviate by 1.5 point and the other 4 by 2.5 points.
//labels = ["12","15","10","13","9.","14","11","16","2","21","4","23","3","24","1","22","20","6.","17","7","18","8","19","5"];
//labels = ["22","3","24","1","21","4","23","2","18","8","19","5","17","7","20","6.","14","11","16","9.","15","10","13","12"];
//labels = ["13","10","15","12","16","11","14","9.","23","4","21","2","22","1","24","3","5","19","8","18","7","17","6.","20"];
//labels = ["6.","20","7","17","5","19","8","18","14","12","15","9.","13","11","16","10","2","24","3","21","4","22","1","23"];
//labels = ["19","5","20","6.","18","8","17","7","12","14","11","13","9.","15","10","16","21","3","24","2","23","1","22","4"];
//labels = ["6.","19","8","17","5","20","7","18","2","24","3","21","1","23","4","22","14","11","16","9.","15","10","13","12"];
//labels = ["21","2","24","3","23","4","22","1","13","11","16","10","14","12","15","9.","18","5","19","8","17","6.","20","7"];

// this labeling has balanced halves; vertices of degree 4 are balanced, 4 vertices of degree 3 deviate by .5 (minimal) and the 4 others by 1.5; opposed faces deviate by 1 (twice), 3 (six times) or 5 (four times) points.
//labels = ["15","9.","16","10","14","12","13","11","20","1","24","5","23","6.","19","2","7","17","4","22","3","21","8","18"];

//this labeling has balanced opposed faces; two halves are balanced, two deviate by 2 points, two deviate by 4 points; 14 vertices of degree 4 are balanced, two deviate by 2 points and two by 4 points; vertices of degree 3 deviate by 0.5, 1.5, 2.5 or 3.5 (two each)
//labels = ["13","10","14","9.","16","11","15","12","22","4","24","2","23","1","21","3","5","18","8","19","7","20","6.","17"];
//labels = ["7","17","5","19","6.","20","8","18","1","22","2","21","4","23","3","24","15","9.","14","12","16","10","13","11"];

// you could also go for this one: opposed faces are balanced; 8 degree 4 vertices are 2 points off, the other 6 are balanced; 6 degree 3 vertices are 2.5 points off and the other 2 are .5 points off; halves are all 2 points off.
//labels = ["18","6.","20","8","17","5","19","7","23","1","21","3","22","4","24","2","12","14","10","16","11","13","9.","15"];

// if you don't care about opposed faces, here is the labeling for you: halves are balanced, vertices of degree 4 are balanced, vertices of degree deviate by .5 except for two which deviate by 1.5.
labels = ["20","5","12","13","19","6.","11","14","10","1","16","23","24","15","2","9.","22","3","8","17","7","18","21","4"];

// coordinate source:  http://dmccooey.com/polyhedra/DeltoidalIcositetrahedron.txt

vertices = [
[ 0.0,  0.0,   C1],
[ 0.0,  0.0,  -C1],
[  C1,  0.0,  0.0],
[ -C1,  0.0,  0.0],
[ 0.0,   C1,  0.0],
[ 0.0,  -C1,  0.0],
[ 1.0,  0.0,  1.0],
[ 1.0,  0.0, -1.0],
[-1.0,  0.0,  1.0],
[-1.0,  0.0, -1.0],
[ 1.0,  1.0,  0.0],
[ 1.0, -1.0,  0.0],
[-1.0,  1.0,  0.0],
[-1.0, -1.0,  0.0],
[ 0.0,  1.0,  1.0],
[ 0.0,  1.0, -1.0],
[ 0.0, -1.0,  1.0],
[ 0.0, -1.0, -1.0],
[  C0,   C0,   C0],
[  C0,   C0,  -C0],
[  C0,  -C0,   C0],
[  C0,  -C0,  -C0],
[ -C0,   C0,   C0],
[ -C0,   C0,  -C0],
[ -C0,  -C0,   C0],
[ -C0,  -C0,  -C0]];
faces = [
[ 14 , 18,  6,  0],
[  8 , 22, 14,  0],
[ 16 , 24,  8,  0],
[  6 , 20, 16,  0],
[ 15 , 23,  9,  1],
[  7 , 19, 15,  1],
[ 17 , 21,  7,  1],
[  9 , 25, 17,  1],
[ 10 , 19,  7,  2],
[  6 , 18, 10,  2],
[ 11 , 20,  6,  2],
[  7 , 21, 11,  2],
[ 12 , 22,  8,  3],
[  9 , 23, 12,  3],
[ 13 , 25,  9,  3],
[  8 , 24, 13,  3],
[ 15 , 19, 10,  4],
[ 12 , 23, 15,  4],
[ 14 , 22, 12,  4],
[ 10 , 18, 14,  4],
[ 16 , 20, 11,  5],
[ 13 , 24, 16,  5],
[ 17 , 25, 13,  5],
[ 11 , 21, 17,  5]];

function add3(v, i = 0, r) = 
    i < len(v) ? 
        i == 0 ?
            add3(v, 1, v[0]) :
            add3(v, i + 1, r + v[i]) :
        r;

function facecoord(n) = [for(i = [0 : len(faces[n]) - 1]) vertices[faces[n][i]] ]; // returns list of the coordinates of the vertices on this face

module facetext(vert,txt) {
    bar = add3(vert)/len(vert); // barycentre
    length = norm(bar);     // radial distance
    b = acos(bar.z/length); // inclination angle
    c = atan2(bar.y,bar.x); // azimuthal angle
    rotate([0,b,c]) translate([0,0,length +minko])linear_extrude(txt_depth,center=true) text(text=txt,size = txt_size,font=txt_font, halign = "center", valign = "center");
}

scale(scafa)
difference() {
    intersection() {
        minkowski($fn=minkfn){
            polyhedron(points = vertices, faces=faces,  convexity = 20);
            sphere(minko);
        };
        sphere(ordi-roll,$fn=minkfn);
    }
    for(i=[0:len(faces)-1]) facetext(facecoord(i),labels[i]);
}