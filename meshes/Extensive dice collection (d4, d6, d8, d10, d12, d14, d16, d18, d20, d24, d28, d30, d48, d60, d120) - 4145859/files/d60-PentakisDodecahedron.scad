// Pentakis Dodecahedron
// 3 sided faces = 60
// degree 6 vertices = 20
// degree 5 vertices = 12
// halves >= 12 (dodecahedron)

txt_depth = .12;
txt_size = .54;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 32; // scale factor which sets the diameter (distance from one [deg 4]-vertex to the vertex opposite) 
minko = 0.4; //chamfer the edges [0 = disabled]
roll = -0.1; // another effet to round the dice (intersects with a smaller sphere) [0 = disabled]
minkfn = 80; //sets the $fn variable for chamfer and the sphere

C0 = 0.927050983124842272306880251548;
C1 = 1.33058699733550141141687582919;
C2 = 2.15293498667750705708437914596;
C3 = 2.427050983124842272306880251548;

ordi = sqrt(3)*1.5+minko; //original diameter
scafa = diam*.5/ordi; //scale factor

//labels = ["1","2","3","4","5","6.","7","8","9.","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60"]; //default labeling

// opposed faces are balanced and the best balance of the vertices I found: degree 5 vertices have minimal unbalance (.5 points off) while degree 6 vertices are balanced except for four; halves are balanced except for two halves which deviate by 2 points
labels = ["13","32","15","33","60","2","44","30","43","34","59","27","18","31","17","48","1","28","46","29","36","38","19","6.","53","20","26","51","40","16","41","45","21","10","35","25","8","55","42","23","50","58","14","9.","22","37","49","4","56","7","24","54","5","57","12","11","39","52","47","3"]; 

// a good compromise (and also the best balance of the vertices I found): degree 5 vertices have minimal unbalance (.5 points off) while degree 6 vertices are balanced except for four; opposed faces are all balanced except for four pairs which deviate by 1; halves are also balanced except for four which deviate by 1.
//labels = ["13","31","15","33","60","2","43","29","44","34","59","27","17","32","18","48","1","28","46","30","36","38","19","6.","53","20","26","51","39","16","41","45","22","9.","35","25","8","54","42","23","50","58","14","10","21","37","49","4","56","7","24","55","5","57","12","11","40","52","47","3"];

// another good compromise: degree 5 vertices have minimal unbalance (.5 points off) while degree 6 vertices are balanced except for six; opposed faces are all balanced except for two pairs which deviate by 1; halves are also balanced except for four which deviate by 1.
//labels = ["13","31","15","33","60","2","44","30","43","34","59","27","18","32","17","48","1","28","46","29","36","38","19","6.","53","20","26","51","40","16","41","45","21","10","35","25","8","55","42","23","50","58","14","9.","22","37","49","4","56","7","24","54","5","57","12","11","39","52","47","3"];

// halves are balanced; six pairs of opposite faces deviate by 1 point; degree 5 vertices are as balanced as possible; degree 6 vertices are balanced except for six which deviate by 1 point.
//labels = ["13","31","15","33","60","2","43","30","44","34","59","27","17","32","18","48","1","28","46","29","36","38","19","6.","53","20","26","51","39","16","41","45","21","10","35","25","8","54","42","23","50","58","14","9.","22","37","49","4","56","7","24","55","5","57","12","11","40","52","47","3"];


// the 12 halves are inherited from those of the underlying dodecahedron. As such, balancing the vertices of degree 5 also balances the halves. The degree 5 vertex  cannot be perfectly balanced though (there are always at leats .5 off). But the halves could be. There might be other reasonnable ways of defining halves.

// source:  http://dmccooey.com/polyhedra/PentakisDodecahedron.txt

vertices = [
[ 0.0,   C0,   C3],
[ 0.0,   C0,  -C3],
[ 0.0,  -C0,   C3],
[ 0.0,  -C0,  -C3],
[  C3,  0.0,   C0],
[  C3,  0.0,  -C0],
[ -C3,  0.0,   C0],
[ -C3,  0.0,  -C0],
[  C0,   C3,  0.0],
[  C0,  -C3,  0.0],
[ -C0,   C3,  0.0],
[ -C0,  -C3,  0.0],
[  C1,  0.0,   C2],
[  C1,  0.0,  -C2],
[ -C1,  0.0,   C2],
[ -C1,  0.0,  -C2],
[  C2,   C1,  0.0],
[  C2,  -C1,  0.0],
[ -C2,   C1,  0.0],
[ -C2,  -C1,  0.0],
[ 0.0,   C2,   C1],
[ 0.0,   C2,  -C1],
[ 0.0,  -C2,   C1],
[ 0.0,  -C2,  -C1],
[ 1.5,  1.5,  1.5],
[ 1.5,  1.5, -1.5],
[ 1.5, -1.5,  1.5],
[ 1.5, -1.5, -1.5],
[-1.5,  1.5,  1.5],
[-1.5,  1.5, -1.5],
[-1.5, -1.5,  1.5],
[-1.5, -1.5, -1.5]];
faces = [
[  2 ,  0, 12],
[ 26 ,  2, 12],
[  4 , 26, 12],
[ 24 ,  4, 12],
[  0 , 24, 12],
[  1 ,  3, 13],
[ 25 ,  1, 13],
[  5 , 25, 13],
[ 27 ,  5, 13],
[  3 , 27, 13],
[  0 ,  2, 14],
[ 28 ,  0, 14],
[  6 , 28, 14],
[ 30 ,  6, 14],
[  2 , 30, 14],
[  3 ,  1, 15],
[ 31 ,  3, 15],
[  7 , 31, 15],
[ 29 ,  7, 15],
[  1 , 29, 15],
[  5 ,  4, 16],
[ 25 ,  5, 16],
[  8 , 25, 16],
[ 24 ,  8, 16],
[  4 , 24, 16],
[  4 ,  5, 17],
[ 26 ,  4, 17],
[  9 , 26, 17],
[ 27 ,  9, 17],
[  5 , 27, 17],
[  6 ,  7, 18],
[ 28 ,  6, 18],
[ 10 , 28, 18],
[ 29 , 10, 18],
[  7 , 29, 18],
[  7 ,  6, 19],
[ 31 ,  7, 19],
[ 11 , 31, 19],
[ 30 , 11, 19],
[  6 , 30, 19],
[ 10 ,  8, 20],
[ 28 , 10, 20],
[  0 , 28, 20],
[ 24 ,  0, 20],
[  8 , 24, 20],
[  8 , 10, 21],
[ 25 ,  8, 21],
[  1 , 25, 21],
[ 29 ,  1, 21],
[ 10 , 29, 21],
[  9 , 11, 22],
[ 26 ,  9, 22],
[  2 , 26, 22],
[ 30 ,  2, 22],
[ 11 , 30, 22],
[ 11 ,  9, 23],
[ 31 , 11, 23],
[  3 , 31, 23],
[ 27 ,  3, 23],
[  9 , 27, 23]];

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