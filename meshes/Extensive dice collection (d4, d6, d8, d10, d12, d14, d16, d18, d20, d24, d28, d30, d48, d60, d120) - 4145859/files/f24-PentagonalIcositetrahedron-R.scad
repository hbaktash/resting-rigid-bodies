//Pentagonal Icositetrahedron (right-handed)
// 5 sided faces = 24
// degree 4 vertices = 6
// degree 3 vertices = 32
// halves = 8

txt_depth = .25;
txt_size = .45;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 24; // scale factor which sets the diameter (distance from one [deg 4]-vertex to the vertex opposite) 
minko = 0.2; //chamfer the edges [0 = disabled]
roll = 0.; // another effet to round the dice (intersects with a smaller sphere) [0 = disabled]
minkfn = 80; //sets the $fn variable for chamfer and the sphere

C0 = 0.218796643000480441021525033712;
C1 = 0.740183741369857222808508165943;
C2 = 1.02365617811269018236347687527;
C3 = 1.36141015192644253450100912043;

ordi = C3+minko; //original diameter
scafa = diam*.5/ordi; //scale factor

//labels = ["1","2","3","4","5","6.","7","8","9.","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24"]; // default ordering

//The following labelling has balanced halves and vertices of degree 4. Vertices of degree 3 deviate by 10.5 (once), 9.5 (twice), 8.5 (four times), 7.5 (thrice), 6.5 (twice), 5.5 (twice), 4.5 (four times), 3.5 (thrice), 2.5 (twice), 1.5 (four times) and .5 (five times)
labels = ["3","19","21","7","24","8","4","14","1","9.","23","17","12","20","16","2","18","15","11","6.","5","10","13","22"];


////The following labelling has balanced halves; vertices of degree 4 deviate by 9 points (four times) or are balanced (twice); vertices of degree 3 deviate by 7.5 (once), 5.5 (five times), 4.5 (twice), 3.5 (thrice), 2.5 (seven times), 1.5 (five times) and .5 (nine times)
//labels = ["17","8","13","3","16","9.","12","4","7","19","2","22","5","20","1","24","21","6.","18","14","23","11","15","10"];

// coordinate source:  http://dmccooey.com/polyhedra/RpentagonalIcositetrahedron.txt
vertices = [
[0.0, 0.0,  C3],
[0.0, 0.0, -C3],
[ C3, 0.0, 0.0],
[-C3, 0.0, 0.0],
[0.0,  C3, 0.0],
[0.0, -C3, 0.0],
[ C1, -C0,  C2],
[ C1,  C0, -C2],
[-C1,  C0,  C2],
[-C1, -C0, -C2],
[ C2, -C1,  C0],
[ C2,  C1, -C0],
[-C2,  C1,  C0],
[-C2, -C1, -C0],
[ C0, -C2,  C1],
[ C0,  C2, -C1],
[-C0,  C2,  C1],
[-C0, -C2, -C1],
[ C0,  C1,  C2],
[ C0, -C1, -C2],
[-C0, -C1,  C2],
[-C0,  C1, -C2],
[ C2,  C0,  C1],
[ C2, -C0, -C1],
[-C2, -C0,  C1],
[-C2,  C0, -C1],
[ C1,  C2,  C0],
[ C1, -C2, -C0],
[-C1, -C2,  C0],
[-C1,  C2, -C0],
[ C1,  C1,  C1],
[ C1,  C1, -C1],
[ C1, -C1,  C1],
[ C1, -C1, -C1],
[-C1,  C1,  C1],
[-C1,  C1, -C1],
[-C1, -C1,  C1],
[-C1, -C1, -C1]];
faces = [
[ 18 , 30, 22,  6,  0],
[  8 , 34, 16, 18,  0],
[ 20 , 36, 24,  8,  0],
[  6 , 32, 14, 20,  0],
[ 19 , 33, 23,  7,  1],
[  9 , 37, 17, 19,  1],
[ 21 , 35, 25,  9,  1],
[  7 , 31, 15, 21,  1],
[ 23 , 33, 27, 10,  2],
[ 11 , 31,  7, 23,  2],
[ 22 , 30, 26, 11,  2],
[ 10 , 32,  6, 22,  2],
[ 25 , 35, 29, 12,  3],
[ 13 , 37,  9, 25,  3],
[ 24 , 36, 28, 13,  3],
[ 12 , 34,  8, 24,  3],
[ 29 , 35, 21, 15,  4],
[ 16 , 34, 12, 29,  4],
[ 26 , 30, 18, 16,  4],
[ 15 , 31, 11, 26,  4],
[ 28 , 36, 20, 14,  5],
[ 17 , 37, 13, 28,  5],
[ 27 , 33, 19, 17,  5],
[ 14 , 32, 10, 27,  5]];

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
    rotate([180,b,c]) translate([0,0,length +minko]) linear_extrude(txt_depth,center=true) text(text=txt,size = txt_size, halign = "center", valign = "center");
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