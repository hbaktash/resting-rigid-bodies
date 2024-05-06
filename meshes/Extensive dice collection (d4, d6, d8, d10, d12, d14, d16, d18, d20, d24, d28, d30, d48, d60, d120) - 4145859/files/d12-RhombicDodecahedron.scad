// Rhombic Dodecahedron 
// 4 sided faces = 12
// degree 3 vertices =6
// degree 4 vertices =8

txt_depth = .09;
txt_size = .475; //.55
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 24; // scale factor which sets the diameter (distance from one [deg 4]-vertex to the vertex opposite) 
minko = 0; //.4 //chamfer the edges
roll = 0.195; // 0 // another effet to round the dice (intersects with a smaller sphere) [ disabled when <0]
minkfn = 80; //sets the $fn variable for chamfer and the sphere

C0 = 0.530330085889910643300633271579;
C1 = 1.06066017177982128660126654316;

ordi = C1+minko; //original diameter
scafa = diam*.5/ordi; //scale factor

//labels = ["1","2","3","4","5","6.","7","8","9.","10","11","12"]; // default labeling

// these three labeling minimise vertex imbalance, and have balanced opposed faces. (degree 3 vertices cannot be better balanced)
//labels = ["1","11","7","4","10","6.","12","3","5","9.","2","8"];
//labels = ["1","11","7","5","9.","6.","12","4","3","8","2","10"];
//labels = ["1","11","8","6.","9.","5","12","4","3","7","2","10"];

// these two labeling have balanced degree 4 vertices, the imbalance by degree 3 vertices is minimised. The opposed faces are no longer balanced.
//labels = ["1","7","12","10","8","3","11","5","2","4","6.","9."];
//labels = ["1","12","6.","4","9.","7","10","3","5","11","2","8"];

labels = ["1", "11", "8", "6.", "9.", "5", "12", "4", "3", "7", "2", "10"];

// I could not find a way to define halves for this solid

// coordinate source:  http://dmccooey.com/polyhedra/RhombicDodecahedron.txt

vertices = [
[0.0, 0.0,  C1],
[0.0, 0.0, -C1],
[ C1, 0.0, 0.0],
[-C1, 0.0, 0.0],
[0.0,  C1, 0.0],
[0.0, -C1, 0.0],
[ C0,  C0,  C0],
[ C0,  C0, -C0],
[ C0, -C0,  C0],
[ C0, -C0, -C0],
[-C0,  C0,  C0],
[-C0,  C0, -C0],
[-C0, -C0,  C0],
[-C0, -C0, -C0]];
faces = [
[  2 ,  8,  0,  6],
[  4 ,  7,  2,  6],
[  0 , 10,  4,  6],
[  2 ,  7,  1,  9],
[  5 ,  8,  2,  9],
[  1 , 13,  5,  9],
[  3 , 13,  1, 11],
[  4 , 10,  3, 11],
[  1 ,  7,  4, 11],
[  3 , 10,  0, 12],
[  5 , 13,  3, 12],
[  0 ,  8,  5, 12]];




function add3(v, i = 0, r) = 
    i < len(v) ? 
        i == 0 ?
            add3(v, 1, v[0]) :
            add3(v, i + 1, r + v[i]) :
        r;

function facecoord(n) = [for(i = [0 : len(faces[n]) - 1]) vertices[faces[n][i]] ];

module facetext(vert,txt) {
    bar = add3(vert)/len(vert); // barycentre
    length = norm(bar);     // radial distance
    b = acos(bar.z/length); // inclination angle
    c = atan2(bar.y,bar.x); // azimuthal angle
    rotate([0,b,c]) translate([0,0,length +minko])linear_extrude(txt_depth,center=true) text(text=txt,size = txt_size,font = txt_font, halign = "center", valign = "center");
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