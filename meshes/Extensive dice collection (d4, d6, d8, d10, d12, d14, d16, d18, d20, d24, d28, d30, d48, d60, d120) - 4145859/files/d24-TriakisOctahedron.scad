//Triakis Octahedron
// 3 sided faces = 24
// degree 8 vertices = 6
// degree 3 vertices = 8
// halves = 6

txt_depth = .11;
txt_size = .48;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 30; // scale factor which sets the diameter (distance from one [deg 4]-vertex to the vertex opposite) 
minko = 0.2; //chamfer the edges [0 = disabled]
roll = 0.6; // another effet to round the dice (intersects with a smaller sphere) [0 = disabled]
minkfn = 80; //sets the $fn variable for chamfer and the sphere

// source:  http://dmccooey.com/polyhedra/TriakisOctahedron.txt

C0 = 2.41421356237309504880168872421;

ordi = C0 +minko; //original diameter
scafa = diam*.5/ordi; //scale factor

//labels = ["1","2","3","4","5","6.","7","8","9.","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24"]; //default labeling

// There are many labelings which have balanced opposed faces, balanced halves, balanced vertices of degree 8 and minimal imbalance at vertices of degree 3 (that is there are .5 points off). This one tries to be so, that faces with a long edge in common add up to something close to 25 (but they still deviate by 13, 11, 7, 6, 5 and 2 points [twice each])
labels =["15","6.","16","5","21","12","18","17","3","24","2","11","14","23","1","22","8","7","13","4","20","9.","19","10"];

// If you want that faces which share a long edge add number to a number closer to 25, the following labeling is better (sums remain between 20 and 30). The price is that 4 of the vertices of degree 8 are no longer balanced (they are 1 point off).
//["14","1","22","7","19","12","5","23","10","16","4","17","8","21","9.","15","2","20","13","6.","18","3","24","11"];

// Lastly, if you care mostly about these long edges and opposed faces, the following labeling has: balanced opposed faces; the long edges deviate by 1 point (cannot do better, since the only other nuymber which adds to 25 is on the opposed face); halves and vertices of degree 8 deviate by 2 points; vertices of degree 3 deviate by 2.5 points (except two which deviate by .5 points).
//["16","1","18","5","23","12","6.","22","10","14","4","17","8","21","11","15","3","19","13","2","20","7","24","9."];

vertices = [
[ 0.0,  0.0,   C0],
[ 0.0,  0.0,  -C0],
[  C0,  0.0,  0.0],
[ -C0,  0.0,  0.0],
[ 0.0,   C0,  0.0],
[ 0.0,  -C0,  0.0],
[ 1.0,  1.0,  1.0],
[ 1.0,  1.0, -1.0],
[ 1.0, -1.0,  1.0],
[ 1.0, -1.0, -1.0],
[-1.0,  1.0,  1.0],
[-1.0,  1.0, -1.0],
[-1.0, -1.0,  1.0],
[-1.0, -1.0, -1.0]];
faces = [
[  2 ,  0,  6],
[  4 ,  2,  6],
[  0 ,  4,  6],
[  4 ,  1,  7],
[  2 ,  4,  7],
[  1 ,  2,  7],
[  5 ,  0,  8],
[  2 ,  5,  8],
[  0 ,  2,  8],
[  2 ,  1,  9],
[  5 ,  2,  9],
[  1 ,  5,  9],
[  4 ,  0, 10],
[  3 ,  4, 10],
[  0 ,  3, 10],
[  3 ,  1, 11],
[  4 ,  3, 11],
[  1 ,  4, 11],
[  3 ,  0, 12],
[  5 ,  3, 12],
[  0 ,  5, 12],
[  5 ,  1, 13],
[  3 ,  5, 13],
[  1 ,  3, 13]];

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
    rotate([0,b,c]) translate([0,0,length +minko])linear_extrude(txt_depth,center=true) text(text=txt,size = txt_size, font = txt_font, halign = "center", valign = "center");
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
