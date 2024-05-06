// Pentagonal Trapezohedron
// 4 sided faces = 8
// degree 3 vertices = 8
// degree 4 vertices = 2
// halves = 8

txt_depth = .25;
txt_size = .55;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 24; // scale factor which sets the diameter (distance from one vertex to the vertex opposite) 
minko = 0.45; //rounds the edges
roll = 0.285; // another effet to round the dice (intersects with a smaller sphere) [0 = disabled]
minkfn = 80; //sets the $fn variable for the "rounding"s

C0 = 0.174155349874503261843187247164;
C1 = 0.7071067811865475244008443621048;
C2 = 1.01505176512821780487431272340;

ordi = C2+minko;
scafa = diam*.5/ordi;

//labels = ["1","2","3","4","5","6","7","8"]; // default labeling

//there are two types of nice labelling. Both have balanced vertices of degree 4 ...
// ... some have vertices of degree 3 which deviate by 1.5 (four times) or .5 (five times)
labels= ["1","4","7","6","5","2","3","8"];
//labels= ["1","4","8","5","6","2","3","7"];

// ... and some have vertices of degree 3 which deviate by 2.5 (once), 1.5 (twice) and .5 (five times)
//labels= ["1","4","6","7","5","2","3","8"];
//labels= ["1","4","8","5","6","3","2","7"];

// labellings above try to minimise imbalance of the "halves". Halves are given by looking at an edge which does not have the summit of a pyramid as an end. These edges touch 4 faces.

//this labeling has balanced vertices of degree 4; vertices of degree 3 deviate by 3.5, 2.5, 1.5 or .5 points (twice each); halves deviate by 1 or 0 points (four times each).
//labels = ["1","8","2","7","6","4","5","3"];

// coordinate source: http://dmccooey.com/polyhedra/TetragonalTrapezohedron.txt

vertices = [
[ 0.0,  0.0,  C2],
[ 0.0,  0.0, -C2],
[  C1,  0.0,  C0],
[ -C1,  0.0,  C0],
[ 0.0,   C1,  C0],
[ 0.0,  -C1,  C0],
[ 0.5,  0.5, -C0],
[ 0.5, -0.5, -C0],
[-0.5,  0.5, -C0],
[-0.5, -0.5, -C0]];

faces = [
[ 4 , 6, 2, 0],
[ 3 , 8, 4, 0],
[ 5 , 9, 3, 0],
[ 2 , 7, 5, 0],
[ 7 , 2, 6, 1],
[ 9 , 5, 7, 1],
[ 8 , 3, 9, 1],
[ 6 , 4, 8, 1]];


function add3(v, i = 0, r) = 
    i < len(v) ? 
        i == 0 ?
            add3(v, 1, v[0]) :
            add3(v, i + 1, r + v[i]) :
        r;

function facecoord(n) = [for(i = [0 : len(faces[n]) - 1]) vertices[faces[n][i]] ];

//echo(add3(facecoord(4))/4);

module facetext(vert,txt,tv) {
    barpo = add3(vert)/len(vert); // barycentre
    bar = barpo;
    length = norm(bar);     // radial distance
    b = acos(bar.z/length); // inclination angle
    c = atan2(bar.y,bar.x); // azimuthal angle
    lenfac = (norm([C1,0,C2-C0])/2+minko)/norm(barpo);
    barpof = [for(j=[0:2]) barpo[j]*lenfac*-1]; // stretching coordiante to compensate for chamfering
    translate(barpof) rotate([180,b,c]) linear_extrude(txt_depth,center=false) text(text=txt,size = txt_size, font=txt_font, halign = "center", valign = "center");
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
    for(i=[0:len(faces)-1]) facetext(facecoord(i),labels[i],faces[i][0]);
}