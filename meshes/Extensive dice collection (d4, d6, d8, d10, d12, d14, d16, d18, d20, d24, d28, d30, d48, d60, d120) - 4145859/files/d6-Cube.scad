// Cube (Hexahedron)
// 4 sided faces = 6
// degree 3 vertices = 8

txt_depth = .09;
txt_size = .35;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 26; // scale factor which sets the diameter (distance from one [deg 4]-vertex to the vertex opposite) 
minko = 0; //chamfer the edges [0 = disabled]
roll = .16; // another effet to round the dice (intersects with a smaller sphere) [ disabled when <0]
minkfn = 80; //sets the $fn variable for chamfer and the sphere

ordi = sqrt(3)*.5+minko; //original diameter
scafa = diam*.5/ordi; //scale factor

//labels = ["1","2","3","4","5","6"]; // default ordering

// This is the usual labeling: opposed faces are balanced; vertices deviate by 4.5, 3.5, 2.5 and 1.5 points (two each).
labels = ["1", "2", "4", "6", "3", "5"];

//This labeling has one pair of opposed faces which are balanced while the other two pairs deviate by 4 points; six vertices are .5 points off and two are 1.5 points off
//labels = ["1","4","6","2","5","3"];

//This labeling has one pair of opposed faces which are balanced while the other two pairs deviate by 3 points; vertices are .5, 1.5, 2.5 and 3.5 points off (two each).
//labels = ["1","2","6","3","4","5"];

//for the cube balaning halves is the same as balancing vertices

// coordinate source:  http://dmccooey.com/polyhedra/Cube.txt

vertices = [
[ 0.5,  0.5,  0.5],
[ 0.5,  0.5, -0.5],
[ 0.5, -0.5,  0.5],
[ 0.5, -0.5, -0.5],
[-0.5,  0.5,  0.5],
[-0.5,  0.5, -0.5],
[-0.5, -0.5,  0.5],
[-0.5, -0.5, -0.5]];
faces = [
[ 4 , 5, 1, 0],
[ 2 , 6, 4, 0],
[ 1 , 3, 2, 0],
[ 6 , 2, 3, 7],
[ 5 , 4, 6, 7],
[ 3 , 1, 5, 7]];

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
    rotate([0,b,c]) translate([0,0,length +minko])linear_extrude(txt_depth,center=true) text(text=txt,size = txt_size,font = txt_font, halign = "center", valign = "center");
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