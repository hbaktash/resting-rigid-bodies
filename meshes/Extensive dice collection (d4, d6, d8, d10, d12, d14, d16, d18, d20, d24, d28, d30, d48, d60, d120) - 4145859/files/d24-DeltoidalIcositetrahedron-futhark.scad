// Deltoidal Icositetrahedron
// 4 sided faces = 24
// degree 4 vertices = 18
// degree 3 vertices = 8
// "halves" = 6

use <./futhark.ttf>

txt_depth = .09;
txt_size = .35;
txt_font = "Futhark:style=regular";
// Warning: a scale factor is applied later (see below)
diam = 24; // scale factor which sets the diameter (distance from one [deg 4]-vertex to the vertex opposite) 
minko = 0.2; //chamfer the edges [0 = disabled]
roll = 0; // another effet to round the dice (intersects with a smaller sphere) [0 = disabled]
minkfn = 80; //sets the $fn variable for chamfer and the sphere

C0 = 0.773459080339013578400241246316;
C1 = 1.41421356237309504880168872421;

ordi = C1+minko; //original diameter
scafa = diam*.5/ordi; //scale factor

labels = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","W","Y","Z"];

//labels =["ᚠ","ᚺ","ᛏ","ᚢ","ᚾ","ᛒ","ᚦ","ᛁ","ᛖ","ᚨ","ᛃ","ᛗ","ᚱ","ᛈ","ᛚ","ᚲ","ᛇ","ᛜ","ᚷ","ᛉ","ᛞ","ᚹ","ᛊ","ᛟ"]; 

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