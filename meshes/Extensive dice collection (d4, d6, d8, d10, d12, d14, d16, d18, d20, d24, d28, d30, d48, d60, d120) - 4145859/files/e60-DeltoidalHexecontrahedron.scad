// Deltoidal Hexecontahedron
// 4 sided faces = 60
// degree 4 vertices = 30
// degree 5 vertices = 12
// degree 3 vertices = 20
// halves = 20 (icosahedron) + 12 (dodecahedron)

txt_depth = .09;
txt_size = .45;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 24; // scale factor which sets the diameter (distance from one [deg 4]-vertex to the vertex opposite) 
minko = 0.4; //chamfer the edges [0 = disabled]
roll = -0.1; // another effet to round the dice (intersects with a smaller sphere) [0 = disabled]
minkfn = 80; //sets the $fn variable for chamfer and the sphere

C0 = 0.690983005625052575897706582817;
C1 = 0.783457635340899531654962439488;
C2 = 1.11803398874989484820458683437;
C3 = 1.20601132958329828273486227812;
C4 = 1.26766108272719625323969951590;
C5 = 1.80901699437494742410229341718;
C6 = 1.95136732208322818153792016770;
C7 = 2.05111871806809578489466195539;
C8 = 2.23606797749978969640917366873;

ordi = sqrt(C3*C3+C6*C6) +minko; //original diameter
scafa = diam*.5/ordi; //scale factor

labels = ["1","2","3","4","5","6.","7","8","9.","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60"]; //default labeling

// coordinate source:  http://dmccooey.com/polyhedra/DeltoidalHexecontahedron.txt

vertices = [
[0.0, 0.0,  C8],
[0.0, 0.0, -C8],
[ C8, 0.0, 0.0],
[-C8, 0.0, 0.0],
[0.0,  C8, 0.0],
[0.0, -C8, 0.0],
[0.0,  C1,  C7],
[0.0,  C1, -C7],
[0.0, -C1,  C7],
[0.0, -C1, -C7],
[ C7, 0.0,  C1],
[ C7, 0.0, -C1],
[-C7, 0.0,  C1],
[-C7, 0.0, -C1],
[ C1,  C7, 0.0],
[ C1, -C7, 0.0],
[-C1,  C7, 0.0],
[-C1, -C7, 0.0],
[ C3, 0.0,  C6],
[ C3, 0.0, -C6],
[-C3, 0.0,  C6],
[-C3, 0.0, -C6],
[ C6,  C3, 0.0],
[ C6, -C3, 0.0],
[-C6,  C3, 0.0],
[-C6, -C3, 0.0],
[0.0,  C6,  C3],
[0.0,  C6, -C3],
[0.0, -C6,  C3],
[0.0, -C6, -C3],
[ C0,  C2,  C5],
[ C0,  C2, -C5],
[ C0, -C2,  C5],
[ C0, -C2, -C5],
[-C0,  C2,  C5],
[-C0,  C2, -C5],
[-C0, -C2,  C5],
[-C0, -C2, -C5],
[ C5,  C0,  C2],
[ C5,  C0, -C2],
[ C5, -C0,  C2],
[ C5, -C0, -C2],
[-C5,  C0,  C2],
[-C5,  C0, -C2],
[-C5, -C0,  C2],
[-C5, -C0, -C2],
[ C2,  C5,  C0],
[ C2,  C5, -C0],
[ C2, -C5,  C0],
[ C2, -C5, -C0],
[-C2,  C5,  C0],
[-C2,  C5, -C0],
[-C2, -C5,  C0],
[-C2, -C5, -C0],
[ C4,  C4,  C4],
[ C4,  C4, -C4],
[ C4, -C4,  C4],
[ C4, -C4, -C4],
[-C4,  C4,  C4],
[-C4,  C4, -C4],
[-C4, -C4,  C4],
[-C4, -C4, -C4]];
faces = [
[ 32 ,  8,  0, 18],
[ 40 , 56, 32, 18],
[ 38 , 10, 40, 18],
[ 30 , 54, 38, 18],
[  0 ,  6, 30, 18],
[ 31 ,  7,  1, 19],
[ 39 , 55, 31, 19],
[ 41 , 11, 39, 19],
[ 33 , 57, 41, 19],
[  1 ,  9, 33, 19],
[ 34 ,  6,  0, 20],
[ 42 , 58, 34, 20],
[ 44 , 12, 42, 20],
[ 36 , 60, 44, 20],
[  0 ,  8, 36, 20],
[ 37 ,  9,  1, 21],
[ 45 , 61, 37, 21],
[ 43 , 13, 45, 21],
[ 35 , 59, 43, 21],
[  1 ,  7, 35, 21],
[ 39 , 11,  2, 22],
[ 47 , 55, 39, 22],
[ 46 , 14, 47, 22],
[ 38 , 54, 46, 22],
[  2 , 10, 38, 22],
[ 40 , 10,  2, 23],
[ 48 , 56, 40, 23],
[ 49 , 15, 48, 23],
[ 41 , 57, 49, 23],
[  2 , 11, 41, 23],
[ 42 , 12,  3, 24],
[ 50 , 58, 42, 24],
[ 51 , 16, 50, 24],
[ 43 , 59, 51, 24],
[  3 , 13, 43, 24],
[ 45 , 13,  3, 25],
[ 53 , 61, 45, 25],
[ 52 , 17, 53, 25],
[ 44 , 60, 52, 25],
[  3 , 12, 44, 25],
[ 50 , 16,  4, 26],
[ 34 , 58, 50, 26],
[ 30 ,  6, 34, 26],
[ 46 , 54, 30, 26],
[  4 , 14, 46, 26],
[ 47 , 14,  4, 27],
[ 31 , 55, 47, 27],
[ 35 ,  7, 31, 27],
[ 51 , 59, 35, 27],
[  4 , 16, 51, 27],
[ 48 , 15,  5, 28],
[ 32 , 56, 48, 28],
[ 36 ,  8, 32, 28],
[ 52 , 60, 36, 28],
[  5 , 17, 52, 28],
[ 53 , 17,  5, 29],
[ 37 , 61, 53, 29],
[ 33 ,  9, 37, 29],
[ 49 , 57, 33, 29],
[  5 , 15, 49, 29]];


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

