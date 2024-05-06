// Icosahedron
// 3 sided faces = 20
// degree 5 vertices = 12
// halves = 20

txt_depth = .09;
txt_size = .35;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 24; // scale factor which sets the diameter (distance from one vertex to the vertex opposite) 
minko = 0.3; //chamfer the edges [0 = disabled]
roll = -0.1; // another effet to round the dice (intersects with a smaller sphere) [0 = disabled]
minkfn = 80; //sets the $fn variable for chamfer and the sphere

C0 = 0.809016994374947424102293417183;

ordi = sqrt(C0*C0+.25)+minko; //original diameter
scafa = diam*.5/ordi; //scale factor

//labels = ["1","2","3","4","5","6.","7","8","9.","10","11","12","13","14","15","16","17","18","19","20"]; // default ordering

// this labeling has minimal vertex imbalance (vertices are .5 points off; one cannot do better); opposed faces are balanced.
labels = ["14","20","13","3","2","19","18","8","1","7","12","10","4","5","6.","9.","11","17","16","15"];


// Here are some alternative labelings with the same property (non-exhaustive list!). These labelings have slightly different performances on the balance of the halves (average deviation, average square deviation or maximal deviation minimised)
//labels = ["1","4","14","18","15","6.","3","7","17","20","11","13","19","16","9.","10","8","2","5","12"]; // maximal deviation minimised
//labels = ["1","8","12","15","16","5","6.","9.","13","20","11","7","19","18","4","10","14","2","3","17"]; // square deviation minimised
//labels = ["17","18","10","6.","2","19","15","11","3","4","14","12","1","5","13","7","9.","20","16","8"]; // average deviation minimised
//labels = ["14", "9.", "11", "15", "3", "18", "6.", "10", "12", "7", "13", "20", "4", "5", "19", "8", "1", "17", "16", "2"]; 
//labels = ["7","17","5","12","11","10","9.","16","4","14","19","3","6.","20","8","2","18","15","1","13"];

// one can define "halves" by looking at one face, its neighbours (3 of them), and the neighbours of the neighnours (6 of them). This give one half for each face (so 20 halves).
// This labeling has 6 vertices which deviate by 1.5 points and 6 which deviate by .5 points; the halves are all balanced except two which deviate by 2 points; the opposed faces are not balanced at all
//labels = ["1","11","14","7","20","19","8","13","12","2","9.","18","5","16","3","10","17","6.","15","4"];

//This labeling has 2 vertices which deviate by 1.5 points and 10 which deviate by .5 points; opposed faces are balanced; the halves deviate less than in the first serie (but still a lot).
//labels = ["15","11","12","14","1","20","7","9.","10","6.","13","19","5","3","17","8","2","16","18","4"];

//This labeling has 4 vertices which deviate by 1.5 points and 10 which deviate by .5 points; opposed faces are balanced; the halves deviate less than in the first serie (but still a lot).
//labels = ["1","4","13","19","16","5","2","8","17","20","11","15","18","14","12","10","6.","3","7","9."];

// These labelings have balance opposed faces; neither the vertices nor the halves are really balanced (up to 9 points of deviation) but that's the best compromise I could find
//labels = ["3","11","7","8","16","5","13","14","10","18","2","20","4","15","9.","19","1","17","6.","12"];
//labels = ["1","16","11","7","18","3","14","10","5","20","12","17","8","6.","19","9.","4","13","15","2"];

// same but the halves are more balanced and the vertices less so
//labels = ["1","12","11","6.","16","5","15","10","9.","20","7","19","3","17","8","14","2","18","4","13"];

// coordinate source:  http://dmccooey.com/polyhedra/Icosahedron.txt

vertices = [
[ 0.5,  0.0,   C0],
[ 0.5,  0.0,  -C0],
[-0.5,  0.0,   C0],
[-0.5,  0.0,  -C0],
[  C0,  0.5,  0.0],
[  C0, -0.5,  0.0],
[ -C0,  0.5,  0.0],
[ -C0, -0.5,  0.0],
[ 0.0,   C0,  0.5],
[ 0.0,   C0, -0.5],
[ 0.0,  -C0,  0.5],
[ 0.0,  -C0, -0.5]];
faces = [
[ 10 ,  2,  0],
[  5 , 10,  0],
[  4 ,  5,  0],
[  8 ,  4,  0],
[  2 ,  8,  0],
[ 11 ,  1,  3],
[  7 , 11,  3],
[  6 ,  7,  3],
[  9 ,  6,  3],
[  1 ,  9,  3],
[  7 ,  6,  2],
[ 10 ,  7,  2],
[ 11 ,  7, 10],
[  5 , 11, 10],
[  1 , 11,  5],
[  4 ,  1,  5],
[  9 ,  1,  4],
[  8 ,  9,  4],
[  6 ,  9,  8],
[  2 ,  6,  8]];


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
