//Triakis Icosahedron
// 3 sided faces = 60
// degree 10 vertives = 12
// degree 3 vertices = 20
//halves >= 20 (icosahedron)

txt_depth = .09;
txt_size = .425;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 34; // scale factor which sets the diameter (distance from one [deg 4]-vertex to the vertex opposite) 
minko = 0.2; //chamfer the edges [0 = disabled]
roll = 0.3; // another effet to round the dice (intersects with a smaller sphere) [0 = disabled]
minkfn = 80; //sets the $fn variable for chamfer and the sphere

C0 = 1.04955317926133973822831518963;
C1 = 1.69821271704535895291326075851;
C2 = 1.80901699437494742410229341718;
C3 = 2.74776589630669869114157594814;
C4 = 2.927050983124842272306880251548;


ordi = sqrt(C2*C2+C4*C4) +minko; //original diameter
scafa = diam*.5/ordi; //scale factor

//labels = ["1","2","3","4","5","6.","7","8","9.","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60"]; //default labeling

//This labeling has balanced opposed faces; vertices of degree 3 deviate by the minimal amount possible (.5 points); vertices of degree 10 deviate by 3, 2 or 1 point (four times each); halves deviate by 1 point, except for two which deviate by 3 points.  
labels = ["1","40","51","55","9.","28","6.","33","52","21","60","10","44","15","32","58","30","4","57","31","3","29","46","17","39","48","5","54","2","35","26","59","7","56","13","22","23","50","18","36","8","47","42","12","37","27","45","20","16","41","34","24","19","49","14","25","53","11","43","38"];


// This labeling has balanced opposed faces; vertices of degree 3 deviate by the minimal .5 points except for four which deviate by 1.5 points; vertices of degree 10 are balanced except for 4 which deviate by 1 point; halves deviate by 5 points (six times), 3 points (eigth times) and 1 point (six times).
//labels = ["1","39","51","55","9.","28","6.","33","52","22","60","10","43","16","32","57","30","5","56","31","4","29","45","18","40","50","3","54","2","35","26","59","7","58","11","21","23","49","20","37","8","47","42","13","36","27","46","17","15","44","34","25","19","48","14","24","53","12","41","38"];

//This labeling has balanced vertices of degree 10, except for four which deviate by 1 point; vertices of degree 3 deviate by their minimal .5 except for one which deviates by 1.5; opposed faces deviate by 4 points (once), 3 (once), 2 (once) or 1 point (five times), the remaining 22 are balanced; halves deviate by 4 (twice), 3 (four times), 2 (four times), 1 (eight times) or 0 (twice) points. 
//labels = ["1","40","51","55","9.","28","6.","33","52","21","60","10","44","16","31","57","29","5","56","32","4","30","45","17","41","49","3","54","2","35","26","59","7","58","11","22","23","50","18","36","8","47","42","12","37","27","46","19","15","43","34","24","20","48","13","25","53","14","39","38"];


//this labeling has balanced halves, except for four which deviate by 1 point; vertices of degree 3 deviate by their minimum (.5 points); verties of degree 10  deviate by 4 (once), 3 (thrice), 2 (twice), 1 (thrice) or 0 points (thrice); opposed faces deviate by 2 (once), 1 (twelve times) or 0 points (seventeen times).
//labels = ["1","40","50","55","9.","28","6.","33","53","20","60","11","45","15","32","56","31","4","58","30","3","29","46","17","38","48","5","54","2","35","26","59","7","57","12","22","23","51","18","36","8","47","42","13","37","27","44","21","16","41","34","24","19","49","14","25","52","10","43","39"];



// coordinate source:  http://dmccooey.com/polyhedra/TriakisIcosahedron.txt

vertices = [
[ C2, 0.0,  C4],
[ C2, 0.0, -C4],
[-C2, 0.0,  C4],
[-C2, 0.0, -C4],
[ C4,  C2, 0.0],
[ C4, -C2, 0.0],
[-C4,  C2, 0.0],
[-C4, -C2, 0.0],
[0.0,  C4,  C2],
[0.0,  C4, -C2],
[0.0, -C4,  C2],
[0.0, -C4, -C2],
[0.0,  C0,  C3],
[0.0,  C0, -C3],
[0.0, -C0,  C3],
[0.0, -C0, -C3],
[ C3, 0.0,  C0],
[ C3, 0.0, -C0],
[-C3, 0.0,  C0],
[-C3, 0.0, -C0],
[ C0,  C3, 0.0],
[ C0, -C3, 0.0],
[-C0,  C3, 0.0],
[-C0, -C3, 0.0],
[ C1,  C1,  C1],
[ C1,  C1, -C1],
[ C1, -C1,  C1],
[ C1, -C1, -C1],
[-C1,  C1,  C1],
[-C1,  C1, -C1],
[-C1, -C1,  C1],
[-C1, -C1, -C1]];
faces = [
[  8 ,  0, 12],
[  2 ,  8, 12],
[  0 ,  2, 12],
[  3 ,  1, 13],
[  9 ,  3, 13],
[  1 ,  9, 13],
[  2 ,  0, 14],
[ 10 ,  2, 14],
[  0 , 10, 14],
[ 11 ,  1, 15],
[  3 , 11, 15],
[  1 ,  3, 15],
[  5 ,  0, 16],
[  4 ,  5, 16],
[  0 ,  4, 16],
[  4 ,  1, 17],
[  5 ,  4, 17],
[  1 ,  5, 17],
[  6 ,  2, 18],
[  7 ,  6, 18],
[  2 ,  7, 18],
[  7 ,  3, 19],
[  6 ,  7, 19],
[  3 ,  6, 19],
[  9 ,  4, 20],
[  8 ,  9, 20],
[  4 ,  8, 20],
[ 10 ,  5, 21],
[ 11 , 10, 21],
[  5 , 11, 21],
[  8 ,  6, 22],
[  9 ,  8, 22],
[  6 ,  9, 22],
[ 11 ,  7, 23],
[ 10 , 11, 23],
[  7 , 10, 23],
[  4 ,  0, 24],
[  8 ,  4, 24],
[  0 ,  8, 24],
[  9 ,  1, 25],
[  4 ,  9, 25],
[  1 ,  4, 25],
[ 10 ,  0, 26],
[  5 , 10, 26],
[  0 ,  5, 26],
[  5 ,  1, 27],
[ 11 ,  5, 27],
[  1 , 11, 27],
[  8 ,  2, 28],
[  6 ,  8, 28],
[  2 ,  6, 28],
[  6 ,  3, 29],
[  9 ,  6, 29],
[  3 ,  9, 29],
[  7 ,  2, 30],
[ 10 ,  7, 30],
[  2 , 10, 30],
[ 11 ,  3, 31],
[  7 , 11, 31],
[  3 ,  7, 31]];


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
    rotate([0,b,c]) translate([0,0,length +minko])linear_extrude(txt_depth,center=true) text(text=txt,size = txt_size, font=txt_font,halign = "center", valign = "center");
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
