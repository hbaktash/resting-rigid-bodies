// Dodecahedron
// 5 sided faces = 12
// degree 3 vertices = 20
// "halves" = 12

txt_depth = .09;
txt_size = .75;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 24; // scale factor which sets the diameter (distance from one vertex to the vertex opposite) 
minko = 0.4; //rounds the edges
roll = 0; // another effet to round the dice (intersects with a smaller sphere) [0 = disabled]
minkfn = 80; //sets the $fn variable for the "rounding"s

C0 = 0.809016994374947424102293417183;
C1 = 1.30901699437494742410229341718;

ordi = sqrt(C1*C1+.25)+minko;
scafa = diam*.5/ordi;

//labels = ["1","2","3","4","5","6.","7","8","9.","10","11","12"]; // default labeling


labels = ["1", "7", "11", "4", "12", "3", "6.", "9.", "8", "5", "10", "2"]; // opposed faces are balanced and vertex imbalance minimised: two are 3.5 points off, six are off by 2.5, six are off by 1.5 and six by .5; two halves are off by 7 points, four by 5, two by 3 and four by 1 point.
// it so turns out that minimising vertex imbalance give that the opposed faces are balanced.
// in some sense this is the best choice (unique up to symmetries)

// if you look at one face and its 5 neighbours, this builds a group of 6 faces which one may consider to be a "half".
//The following labeling looses a bit on the balance of the vertices and the opposed faces but has more balanced halves: four vertices are 3.5 points off, four are off by 2.5, six are off by 1.5 and six by .5; four halves are off by 6 points, two by 5, two by 3 and four are balanced; two pairs of opposed faces are off by 2, two by 1 and two are balanced.
//labels = ["10","1","9.","8","3","6.","11","7","2","12","5","4"];

//The following labeling has seriously unbalanced vertices; two opposed faces are off by 10 points, two by 6 and two by 2; but the halves are as balanced as it gets: only two are off by 2 points, the remaining 10 halves are balanced.
//labels = ["10","6.","1","7","9.","4","5","8","11","12","3","2"];



// coordinate source: http://dmccooey.com/polyhedra/Dodecahedron.txt

vertices = [
[ 0.0,  0.5,   C1],
[ 0.0,  0.5,  -C1],
[ 0.0, -0.5,   C1],
[ 0.0, -0.5,  -C1],
[  C1,  0.0,  0.5],
[  C1,  0.0, -0.5],
[ -C1,  0.0,  0.5],
[ -C1,  0.0, -0.5],
[ 0.5,   C1,  0.0],
[ 0.5,  -C1,  0.0],
[-0.5,   C1,  0.0],
[-0.5,  -C1,  0.0],
[  C0,   C0,   C0],
[  C0,   C0,  -C0],
[  C0,  -C0,   C0],
[  C0,  -C0,  -C0],
[ -C0,   C0,   C0],
[ -C0,   C0,  -C0],
[ -C0,  -C0,   C0],
[ -C0,  -C0,  -C0]];
faces = [
[ 12 ,  4, 14,  2,  0],
[ 16 , 10,  8, 12,  0],
[  2 , 18,  6, 16,  0],
[ 17 , 10, 16,  6,  7],
[ 19 ,  3,  1, 17,  7],
[  6 , 18, 11, 19,  7],
[ 15 ,  3, 19, 11,  9],
[ 14 ,  4,  5, 15,  9],
[ 11 , 18,  2, 14,  9],
[  8 , 10, 17,  1, 13],
[  5 ,  4, 12,  8, 13],
[  1 ,  3, 15,  5, 13]];



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