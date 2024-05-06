// Triakis Tetrahedron
// 3 sided faces = 12


txt_depth = .25;
txt_size = .4;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 27; // scale factor which sets the diameter (distance from one [deg 4]-vertex to the vertex opposite) 
minko = 0.4; //chamfer the edges [0 = disabled]
roll = .5; // another effet to round the dice (intersects with a smaller sphere) [disabled if <= 0]
minkfn = 80; //sets the $fn variable for chamfer and the sphere

C0 = 0.636396103067892771960759925894;
C1 = 1.06066017177982128660126654316;

ordi = sqrt(3)*C1+minko; //original diameter
scafa = diam*.5/ordi; //scale factor


//labels = ["1","2","3","4","5","6.","7","8","9.","10","11","12"]; //default labeling

// the following labeling has balanced vertices (at vertices of degree 3 there is a minimal deviation of .5 points) and the number on the edges are also balaced.
labels = ["1","7","11","2","9.","8","5","3","12","6.","10","4"];

// coordinate source:  http://dmccooey.com/polyhedra/TriakisTetrahedron.txt


vertices = [
[ C1,  C1,  C1],
[ C1, -C1, -C1],
[-C1, -C1,  C1],
[-C1,  C1, -C1],
[ C0, -C0,  C0],
[ C0,  C0, -C0],
[-C0,  C0,  C0],
[-C0, -C0, -C0]];
faces = [
[ 2 , 0, 4],
[ 1 , 2, 4],
[ 0 , 1, 4],
[ 1 , 0, 5],
[ 3 , 1, 5],
[ 0 , 3, 5],
[ 3 , 0, 6],
[ 2 , 3, 6],
[ 0 , 2, 6],
[ 2 , 1, 7],
[ 3 , 2, 7],
[ 1 , 3, 7]];

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
    rotate([180,b,c]) translate([0,0,length +minko]) rotate([0,0,-45]) linear_extrude(txt_depth,center=false) text(text=txt,size = txt_size, font=txt_font, halign = "center", valign = "center");
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
