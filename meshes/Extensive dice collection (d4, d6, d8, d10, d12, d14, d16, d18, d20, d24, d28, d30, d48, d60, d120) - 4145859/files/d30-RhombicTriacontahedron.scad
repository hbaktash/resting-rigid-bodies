// Rhombic Triacontahedron
// 4 sided faces = 30
// degree 3 vertices = 20
// degree 5 vertices = 12
// halves = 40

txt_depth = .11;
txt_size = .54;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 28; // scale factor which sets the diameter (distance from one [deg 4]-vertex to the vertex opposite) 
minko = 0.4; //chamfer the edges [0 = disabled]
roll = -.1; // another effet to round the dice (intersects with a smaller sphere) [disabled if <= 0]
minkfn = 80; //sets the $fn variable for chamfer and the sphere


C0 = 0.559016994374947424102293417183;
C1 = 0.904508497187473712051146708591;
C2 = 1.46352549156242113615344012577;

ordi = sqrt(C1*C1+C2*C2)+minko; //original diameter
scafa = diam*.5/ordi; //scale factor


//labels = ["1","2","3","4","5","6.","7","8","9.","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"]; //default labeling

//These labelings have balanced opposed faces; the vertices are as balanced as can be (deviate by .5 points). The halves are not balanced.
labels = ["20","4","15","8","18","26","5","13","23","16","27","11","7","30","2","29","1","24","25","17","9.","3","19","21","10","12","28","22","14","6."]; 
//labels = ["13","6","21","4","15","30","1","16","27","10","25","18","9","24","2","29","7","22","26","20","14","12","8","28","3","23","19","17","11","5"]; 


// here is how to define halves: if you look at the top (or bottom, back, front, right or left) view you see one face and its 4 neighbours. These can be assimilared to the face of the underlying cube. Since this cube has "halves", so does this solid.
// note that there are more than these 8 halves since one could draw the axis through other faces. (Pick a vertex od degree 5, each neighbouring faces will yield a different set of 8 halves)

// In the following labeling are halves quite balanced (32 deviate by .5 points [which is the minimum] and 8 by 1.5 points). The opposed faces and vertices are not balanced at all.
//labels = ["11","29","10","8","5","24","30","2","7","16","25","3","14","28","18","21","17","20","12","1","15","27","4","26","19","9.","22","23","6.","13"];

// The upcoming labelings have fairly balanced halves (deviate up to 3.5) and opposed faces (off by at most 3 points). Vertices are not balanced at all.
//labels = ["13","3","23","12","4","21","11","27","17","7","29","19","20","1","28","5","30","9.","10","26","16","8","25","15","18","6.","22","14","2","24"];
//labels = ["28","17","10","15","19","5","26","14","16","22","12","3","8","30","13","18","1","24","7","11","2","25","9.","29","4","23","6.","27","20","21"];

//The upcoming labelings have fairly balanced halves (up to 2.5 points off) and try to minimise imbalance of the vertices. Opposed faces are not balanced at all.
//labels = ["28","15","8","3","10","13","14","9.","6.","7","18","27","5","30","16","17","29","4","20","21","1","11","24","25","23","26","12","2","22","19"];
//labels = ["19","12","28","21","23","5","3","26","22","27","13","20","25","1","16","15","2","24","29","6.","17","8","9.","14","11","10","7","18","4","30"];

// coordinate source:  http://dmccooey.com/polyhedra/RhombicTriacontahedron.txt

vertices = [
[ C1, 0.0,  C2],
[ C1, 0.0, -C2],
[-C1, 0.0,  C2],
[-C1, 0.0, -C2],
[ C2,  C1, 0.0],
[ C2, -C1, 0.0],
[-C2,  C1, 0.0],
[-C2, -C1, 0.0],
[0.0,  C2,  C1],
[0.0,  C2, -C1],
[0.0, -C2,  C1],
[0.0, -C2, -C1],
[0.0,  C0,  C2],
[0.0,  C0, -C2],
[0.0, -C0,  C2],
[0.0, -C0, -C2],
[ C2, 0.0,  C0],
[ C2, 0.0, -C0],
[-C2, 0.0,  C0],
[-C2, 0.0, -C0],
[ C0,  C2, 0.0],
[ C0, -C2, 0.0],
[-C0,  C2, 0.0],
[-C0, -C2, 0.0],
[ C1,  C1,  C1],
[ C1,  C1, -C1],
[ C1, -C1,  C1],
[ C1, -C1, -C1],
[-C1,  C1,  C1],
[-C1,  C1, -C1],
[-C1, -C1,  C1],
[-C1, -C1, -C1]];
faces = [
[ 14 ,  2, 12,  0],
[ 26 , 10, 14,  0],
[ 16 ,  5, 26,  0],
[ 25 ,  9, 13,  1],
[ 17 ,  4, 25,  1],
[ 27 ,  5, 17,  1],
[ 18 ,  6, 28,  2],
[ 30 ,  7, 18,  2],
[ 14 , 10, 30,  2],
[ 29 ,  6, 19,  3],
[ 13 ,  9, 29,  3],
[ 15 ,  1, 13,  3],
[ 24 ,  8, 20,  4],
[ 16 ,  0, 24,  4],
[ 17 ,  5, 16,  4],
[ 19 ,  6, 18,  7],
[ 31 ,  3, 19,  7],
[ 23 , 11, 31,  7],
[ 28 ,  6, 22,  8],
[ 12 ,  2, 28,  8],
[ 24 ,  0, 12,  8],
[ 22 ,  6, 29,  9],
[ 20 ,  8, 22,  9],
[ 25 ,  4, 20,  9],
[ 23 ,  7, 30, 10],
[ 21 , 11, 23, 10],
[ 26 ,  5, 21, 10],
[ 15 ,  3, 31, 11],
[ 27 ,  1, 15, 11],
[ 21 ,  5, 27, 11]];

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