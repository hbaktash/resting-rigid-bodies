// Hexagonal Dipyramid
// 3 sided faces = 12
// degree 4 vertices = 6
// degree 6 vertices = 2
// halves = 8 (two are the degree 6 vertices)
 
txt_depth = .11;
txt_size = .55;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 27; // scale factor which sets the diameter (distance from one [deg 4]-vertex to the vertex opposite) 
minko = 0.4; //chamfer the edges [0 = disabled]
roll = 0.2; // another effet to round the dice (intersects with a smaller sphere) [disabled if <= 0]
minkfn = 80; //sets the $fn variable for chamfer and the sphere

C0 = 1/sqrt(3);
C1 = 2/sqrt(3);
C2 = C1; // the height of each pyramids in the dipyramid; different from the original coordinates (2)

ordi = max(C2,C1)+minko; //original diameter
scafa = diam*.5/ordi; //scale factor

zint = -C2/3 + 2/(3*C2);//1/sqrt(27); // z intersect of the ray perpendicular to the faces

//labels = ["1","2","3","4","5","6.","7","8","9.","10","11","12"]; //default labeling

//The following labeling has balanced opposed faces; vertices of degree 4 are also balanced; vertices of degree 6 deviate by 1 point; halves deviate by 2 points.
labels = ["1","8","9.","3","6.","11","4","5","12","2","7","10"];

//The following labeling has balanced opposed faces; vertices of degree 4 are also balanced; vertices of degree 6 deviate by 2 points; halves deviate by 1 point.
//labels = ["11","4","5","12","3","6.","8","9.","2","7","10","1"];

//The following labeling has balanced vertices; opposed faces deviate by 1 point; halves deviate by 2 points.
//labels = ["7","9.","4","6.","12","1","10","3","5","11","2","8"];

//The following labeling has balanced vertices and halves; opposed faces deviate by either 2 (four pairs) or 1 point (two pairs).
//labels = ["1","11","7","3","9.","8","5","4","10","6.","2","12"];

//The following labeling has balanced halves although the vertices of degree 6 (which are also halves) deviate by 2 points; the other vertices are balanced; the opposed faces deviate by 1 point.
//labels = ["4","11","6.","3","12","5","8","1","10","7","2","9."];

// coordinate source:  http://dmccooey.com/polyhedra/HexagonalDipyramid.txt

vertices = [
[ 0.0, 0.0,  C2],
[ 0.0, 0.0,  -C2],
[ 1.0,  C0,  0.0],
[ 1.0, -C0,  0.0],
[-1.0,  C0,  0.0],
[-1.0, -C0,  0.0],
[ 0.0,  C1,  0.0],
[ 0.0, -C1,  0.0]];
faces = [
[ 6 , 2, 0],
[ 4 , 6, 0],
[ 5 , 4, 0],
[ 7 , 5, 0],
[ 3 , 7, 0],
[ 2 , 3, 0],
[ 3 , 2, 1],
[ 7 , 3, 1],
[ 5 , 7, 1],
[ 4 , 5, 1],
[ 6 , 4, 1],
[ 2 , 6, 1]];



function add3(v, i = 0, r) = 
    i < len(v) ? 
        i == 0 ?
            add3(v, 1, v[0]) :
            add3(v, i + 1, r + v[i]) :
        r;

function facecoord(n) = [for(i = [0 : len(faces[n]) - 1]) vertices[faces[n][i]] ]; // returns list of the coordinates of the vertices on this face

module facetext(vert,txt,tv) {
    barpo = add3(vert)/len(vert); // barycentre
    bar = add3([barpo,[0,0,( tv == 0 ? 1 : -1)*zint]]);
    length = norm(bar);     // radial distance
    b = acos(bar.z/length); // inclination angle
    c = atan2(bar.y,bar.x); // azimuthal angle
    lenfac = (norm(barpo)+minko)/norm(barpo);
    barpof = [for(j=[0:2]) barpo[j]*lenfac]; // stretching coordiante to compensate for chamfering
    translate(barpof) rotate([0,b,c]) linear_extrude(txt_depth,center=true) text(text=txt,size = txt_size, font=txt_font, halign = "center", valign = "center");
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
    for(i=[0:len(faces)-1]) facetext(facecoord(i),labels[i],faces[i][2]);
}

