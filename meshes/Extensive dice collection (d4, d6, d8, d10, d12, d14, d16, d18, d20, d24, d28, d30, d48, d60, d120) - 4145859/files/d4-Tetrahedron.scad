//Tetrahedron
// 3 sided faces = 4
// degree 3 vertices = 4

txt_depth = .5;
txt_size = .3;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 27; // scale factor which sets the diameter (twice the distance from one vertex to the center) 
minko = 0.2; //chamfer the edges [0 = disabled]
roll = .22; // another effet to round the dice (intersects with a smaller sphere) [ disabled when <0]
minkfn = 80; //sets the $fn variable for chamfer and the sphere

C0 = 0.353553390593273762200422181052;

ordi = sqrt(3)*C0+minko; //original diameter
scafa = diam*.5/ordi; //scale factor


labels = ["1","2","3","4"]; //default labeling

// coordinate source:  http://dmccooey.com/polyhedra/Tetrahedron.txt

vertices = [
[ C0, -C0,  C0],
[ C0,  C0, -C0],
[-C0,  C0,  C0],
[-C0, -C0, -C0]];
faces = [
[ 2 , 1, 0],
[ 3 , 0, 1],
[ 0 , 3, 2],
[ 1 , 2, 3]];


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
    rotate([180,b,c]) translate([0,0,ordi])linear_extrude(txt_depth,center=true) text(text=txt,size = txt_size,font=txt_font, halign = "center", valign = "center");
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