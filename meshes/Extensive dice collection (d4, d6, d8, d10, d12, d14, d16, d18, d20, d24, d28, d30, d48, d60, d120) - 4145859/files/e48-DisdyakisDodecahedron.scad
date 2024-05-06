//Disdyakis Dodecahedron
// a.k.a. Hexakis Oktaeder
// 3 sided faces = 48
// degree 8 vertices = 6
// degree 6 vertices = 8
// degree 4 vertices = 12
// halves = 18

txt_depth = .11;
txt_size = .45;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 32; // scale factor which sets the diameter (distance from one [deg 4]-vertex to the vertex opposite) 
minko = 0.2; //chamfer the edges [0 = disabled]
roll = 0.3; // another effet to round the dice (intersects with a smaller sphere) [0 = disabled]
minkfn = 80; //sets the $fn variable for chamfer and the sphere

C0 = 1.41421356237309504880168872421;
C1 = 1.64075448203408147040144747789;
C2 = 2.67541743733683649131645693113;

ordi = C2+minko; //original diameter
scafa = diam*.5/ordi; //scale factor

labels = ["1","2","3","4","5","6.","7","8","9.","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48"]; //default labeling

// For each axis there are 2 ways to cut the solid in two going only through thus axis (accounting for 6 ways). Then, for each pair of axes, there is another cut. So a total of 9 cuts, yielding 18 possible halves.
// there are also "halves" coming from the underlying cube. But these are automtically balanced when the degree 8 vertices are.

// coordinate source:  http://dmccooey.com/polyhedra/DisdyakisDodecahedron.txt

vertices = [
[0.0, 0.0,  C2],
[0.0, 0.0, -C2],
[ C2, 0.0, 0.0],
[-C2, 0.0, 0.0],
[0.0,  C2, 0.0],
[0.0, -C2, 0.0],
[ C1, 0.0,  C1],
[ C1, 0.0, -C1],
[-C1, 0.0,  C1],
[-C1, 0.0, -C1],
[ C1,  C1, 0.0],
[ C1, -C1, 0.0],
[-C1,  C1, 0.0],
[-C1, -C1, 0.0],
[0.0,  C1,  C1],
[0.0,  C1, -C1],
[0.0, -C1,  C1],
[0.0, -C1, -C1],
[ C0,  C0,  C0],
[ C0,  C0, -C0],
[ C0, -C0,  C0],
[ C0, -C0, -C0],
[-C0,  C0,  C0],
[-C0,  C0, -C0],
[-C0, -C0,  C0],
[-C0, -C0, -C0]];
faces = [
[ 18 ,  6,  0],
[ 14 , 18,  0],
[ 22 , 14,  0],
[  8 , 22,  0],
[ 24 ,  8,  0],
[ 16 , 24,  0],
[ 20 , 16,  0],
[  6 , 20,  0],
[ 21 ,  7,  1],
[ 17 , 21,  1],
[ 25 , 17,  1],
[  9 , 25,  1],
[ 23 ,  9,  1],
[ 15 , 23,  1],
[ 19 , 15,  1],
[  7 , 19,  1],
[ 20 ,  6,  2],
[ 11 , 20,  2],
[ 21 , 11,  2],
[  7 , 21,  2],
[ 19 ,  7,  2],
[ 10 , 19,  2],
[ 18 , 10,  2],
[  6 , 18,  2],
[ 22 ,  8,  3],
[ 12 , 22,  3],
[ 23 , 12,  3],
[  9 , 23,  3],
[ 25 ,  9,  3],
[ 13 , 25,  3],
[ 24 , 13,  3],
[  8 , 24,  3],
[ 19 , 10,  4],
[ 15 , 19,  4],
[ 23 , 15,  4],
[ 12 , 23,  4],
[ 22 , 12,  4],
[ 14 , 22,  4],
[ 18 , 14,  4],
[ 10 , 18,  4],
[ 20 , 11,  5],
[ 16 , 20,  5],
[ 24 , 16,  5],
[ 13 , 24,  5],
[ 25 , 13,  5],
[ 17 , 25,  5],
[ 21 , 17,  5],
[ 11 , 21,  5]];


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
    rotate([0,b,c]) translate([0,0,length +minko])linear_extrude(txt_depth,center=true) text(text=txt,size = txt_size, font=txt_font, halign = "center", valign = "center");
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