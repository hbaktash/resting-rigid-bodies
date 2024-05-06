// Octahedron
// 3 sided faces = 8
// degree 4 vertices = 6

txt_depth = .09;
txt_size = .3;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 24; // scale factor which sets the diameter (distance from one vertex to the vertex opposite) 
minko = 0; //chamfers the edges
roll = 0.2075; // another effet to round the dice (intersects with a smaller sphere) [0 = disabled]
minkfn = 80; //sets the $fn variable for chamfer and the sphere

C0 = 0.7071067811865475244008443621048;

ordi = C0+minko;
scafa = diam*.5/ordi;


//labels = ["1","2","3","4","5","6.","7","8"]; //default labeling
// labels = ["1", "7", "4", "6", "2", "8", "3", "5"]; // opposed faces balanced, vertex imbalance minimised
labels = ["1", "6", "3", "8", "2", "5", "4", "7"] ; // vertices balanced, opposed faces imbalance minimised

// for the Octahedron, balancing halves is the same thing as balancing vertices

// coordinates source:  http://dmccooey.com/polyhedra/Octahedron.txt

vertices = [
[0.0, 0.0,  C0],
[0.0, 0.0, -C0],
[ C0, 0.0, 0.0],
[-C0, 0.0, 0.0],
[0.0,  C0, 0.0],
[0.0, -C0, 0.0]];
faces = [
[ 4 , 2, 0],
[ 3 , 4, 0],
[ 5 , 3, 0],
[ 2 , 5, 0],
[ 5 , 2, 1],
[ 3 , 5, 1],
[ 4 , 3, 1],
[ 2 , 4, 1]];



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
    rotate([0,b,c]) translate([0,0,length + minko])linear_extrude(txt_depth,center=true) text(text=txt,size = txt_size,font=txt_font, halign = "center", valign = "center");
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
