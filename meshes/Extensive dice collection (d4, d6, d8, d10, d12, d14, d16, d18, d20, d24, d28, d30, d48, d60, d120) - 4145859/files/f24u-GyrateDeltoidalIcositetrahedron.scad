//"Gyrate Deltoidal Icositetrahedron",
// 4 sided faces = 24
// degree 4 vertices = 18
// degree 3 vertices = 8
// halves = 2 (plus 8 inherited from cube)


txt_depth = .24;
txt_size = .4;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 24; // scale factor which sets the diameter (distance from one [deg 4]-vertex to the vertex opposite) 
minko = 0.4; //chamfer the edges [0 = disabled]
roll = -0.1; // another effet to round the dice (intersects with a smaller sphere) [0 = disabled]
minkfn = 80; //sets the $fn variable for chamfer and the sphere

C0 = 0.7071067811865475244008443621048;
C1 = 0.773459080339013578400241246316;
C2 = 1.09383632135605431360096498526;
C3 = 1.41421356237309504880168872421;

ordi = sqrt(C2*C2+C1*C1) +minko; //original diameter
scafa = diam*.5/ordi; //scale factor

//labels = ["1","2","3","4","5","6.","7","8","9.","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24"]; //default labeling

// Note that some face are opposed, but not all. So balance of opposed face refers only to those face which make sense

//This labeling has balanced opposed faces and halves; vertices of degree 3 deviate by 1.5 or .5 points (4 times each); vertices of degree 4 deviate by 4 (four times), 3 (once), 2 (four times), 1 (once) and 0 (eight times) points.
labels = ["14","12","13","11","9.","15","10","16","4","22","1","23","3","21","2","24","6.","20","7","17","5","19","8","18"];

//For this labeling opposed faces were ignored (so they deviate a lot). Halves are balanced; vertices of degree 4 deviate by 2 (once), 1 (six times) or 0 (eleven times); vertices of degree 3 deviate by 2.5 (twice), 1.5 (once) and .5 (five times).
//labels = ["16","19","6.","9.","20","12","5","13","23","2","17","8","14","10","22","4","18","7","24","1","21","3","15","11"];


//This labeling has the best vertex balance I could find. Opposed faes are not balanced at all; vertices of degree 4 deviate by 1 (eight times) or 0 (ten times) points; vertices of degree 3 deviate by 2.5 (once), 1.5 (twice) or .5 (five times); halves deviate by 1 (six times) or 0 (four times).
//labels = ["16","19","6.","9.","20","12","5","13","24","1","17","7","14","10","22","4","18","8","23","2","21","3","15","11"];




vertices = [
[ 0.0,  0.0,   C3],
[ 0.0,  0.0,  -C3],
[  C3,  0.0,  0.0],
[ -C3,  0.0,  0.0],
[ 0.0,   C3,  0.0],
[ 0.0,  -C3,  0.0],
[ 0.0,   C2,   C1],
[ 0.0,  -C2,   C1],
[  C2,  0.0,   C1],
[ -C2,  0.0,   C1],
[ 1.0,  0.0, -1.0],
[-1.0,  0.0, -1.0],
[ 1.0,  1.0,  0.0],
[ 1.0, -1.0,  0.0],
[-1.0,  1.0,  0.0],
[-1.0, -1.0,  0.0],
[ 0.0,  1.0, -1.0],
[ 0.0, -1.0, -1.0],
[  C0,   C0,  1.0],
[  C0,  -C0,  1.0],
[ -C0,   C0,  1.0],
[ -C0,  -C0,  1.0],
[  C1,   C1,  -C1],
[  C1,  -C1,  -C1],
[ -C1,   C1,  -C1],
[ -C1,  -C1,  -C1]];
faces = [
[ 20 ,  6, 18,  0],
[ 21 ,  9, 20,  0],
[ 19 ,  7, 21,  0],
[ 18 ,  8, 19,  0],
[ 17 , 23, 10,  1],
[ 11 , 25, 17,  1],
[ 16 , 24, 11,  1],
[ 10 , 22, 16,  1],
[ 13 , 19,  8,  2],
[ 10 , 23, 13,  2],
[ 12 , 22, 10,  2],
[  8 , 18, 12,  2],
[ 14 , 20,  9,  3],
[ 11 , 24, 14,  3],
[ 15 , 25, 11,  3],
[  9 , 21, 15,  3],
[ 12 , 18,  6,  4],
[ 16 , 22, 12,  4],
[ 14 , 24, 16,  4],
[  6 , 20, 14,  4],
[ 15 , 21,  7,  5],
[ 17 , 25, 15,  5],
[ 13 , 23, 17,  5],
[  7 , 19, 13,  5]];


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
    rotate([180,b,c]) translate([0,0,length-.05 +minko]) linear_extrude(txt_depth,center=false) text(text=txt,size = txt_size, font=txt_font,halign = "center", valign = "center");
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
