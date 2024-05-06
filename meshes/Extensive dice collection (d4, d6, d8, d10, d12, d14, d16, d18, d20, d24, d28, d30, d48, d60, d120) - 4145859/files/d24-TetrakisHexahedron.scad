// TetrakisHexahedron
// 3 sided faces = 24
// degree 4 vertices = 6
// degree 6 vertices = 8
// halves = 8

txt_depth = .09;
txt_size = .5;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 26; // scale factor which sets the diameter (distance from one [deg 6]-vertex to the vertex opposite) 
minko = 0.4; //chamfer the edges [0 = disabled]
roll = 0.225; // another effet to round the dice (intersects with a smaller sphere) [0 = disabled]
minkfn = 80; //sets the $fn variable for chamfer and the sphere

C0 = 1.06066017177982128660126654316;
C1 = 1.59099025766973192990189981474;


ordi = sqrt(3)*C0+minko; //original diameter
scafa = diam*.5/ordi; //scale factor

//labels = ["1","2","3","4","5","6.","7","8","9.","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24"]; // default labeling

// this solid can be cut in halves (there are 6 possibilities each yielding 2 halves), and so one try to balance these halves. 
// When the opposed faces are balanced, one half and its opposed half always deviate from the same amount (with opposite sign) from the optimal value 

// In the upcoming 3 labelings vertices and opposed faces are balanced, halves deviate 12, 10, 6 [two times] 4 and 2 points [from 150]); 
labels = ["3","19","11","17","6.","22","8","14","9.","10","13","18","15","16","7","12","4","20","2","24","5","21","1","23"]; //
//labels = ["5","23","1","21","2","20","4","24","18","13","10","9.","12","7","16","15","8","22","6.","14","3","17","11","19"];
//labels = ["14","6.","22","8","19","11","17","3","13","10","9.","18","15","12","7","16","2","24","4","20","1","23","5","21"]; // similar to (3-19-11-17...) above

// some labelings are very similar to others (but cannot be obtained by just turning the dice around or mirror symmetry)

// this one is a good compromise: all vertices but 4 are balanced (those which are not balanced are 1 point off), 6 from the 12 halves are balanced (the others are 4 points off) and the opposed faces are all 1 point off.
//labels = ["22","1","24","3","23","4","21","2","17","6.","19","8","20","7","18","5","10","14","11","15","12","16","9.","13"];

//These are my personal favourites: vertices and halves are balanced, the opposed vertices are either 1 or 3 points off
//labels = ["16","9.","14","11","13","12","15","10","5","18","7","20","6.","17","8","19","22","3","24","1","23","2","21","4"];
//labels = ["11","14","12","13","10","15","9.","16","18","8","17","7","20","6.","19","5","2","23","1","24","3","22","4","21"];

//These have balanced halves and opposed faces; 2 vertices are 3 points off, 6 are 2 points off and 6 others are 1 point off.
//labels = ["24","4","22","2","21","1","23","3","18","6.","20","8","19","7","17","5","13","11","15","9.","14","12","16","10"];
//labels = ["1","23","3","21","2","24","4","22","12","16","10","14","9.","13","11","15","18","6.","20","8","19","7","17","5"]
//labels = ["12","14","10","16","11","13","9.","15","1","21","3","23","4","24","2","22","19","7","17","5","18","6.","20","8"]

// note that one could define further halves by looking at three adjacent faces of the underlying cube. However if the degree 4 vertices are balanced, then these halves are automatically balanced.

// coordinate source:  http://dmccooey.com/polyhedra/TetrakisHexahedron.txt


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
[ 10 ,  6,  0],
[ 12 , 10,  0],
[  8 , 12,  0],
[  6 ,  8,  0],
[  9 ,  7,  1],
[ 13 ,  9,  1],
[ 11 , 13,  1],
[  7 , 11,  1],
[  8 ,  6,  2],
[  9 ,  8,  2],
[  7 ,  9,  2],
[  6 ,  7,  2],
[ 11 , 10,  3],
[ 13 , 11,  3],
[ 12 , 13,  3],
[ 10 , 12,  3],
[  7 ,  6,  4],
[ 11 ,  7,  4],
[ 10 , 11,  4],
[  6 , 10,  4],
[ 12 ,  8,  5],
[ 13 , 12,  5],
[  9 , 13,  5],
[  8 ,  9,  5]];

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
    rotate([0,b,c]) translate([0,0,length +minko])linear_extrude(txt_depth,center=true) text(text=txt,size = txt_size, font =txt_font, halign = "center", valign = "center");
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
