// Octahedron
// 4 sided faces = 8
// degree 4 vertices = 2
// degree 3 vertices = 8

txt_depth = .09;
txt_size = .3;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 24; // scale factor which sets the diameter (distance from one vertex to the vertex opposite) 
minko = 0.1; //chamfers the edges
roll = -0.3; // another effet to round the dice (intersects with a smaller sphere) [0 = disabled]
minkfn = 80; //sets the $fn variable for chamfer and the sphere

C0 = 1;
C1 = 0.6;
C2 = C1*(1-C1)/(1+C1);


ordi = C0+minko;
scafa = diam*.5/ordi;


labels = ["1","2","3","4","5","6.","7","8","9","10"]; //default labeling
// labels = ["1", "7", "4", "6", "2", "8", "3", "5"]; // opposed faces balanced, vertex imbalance minimised
//labels = ["1", "6", "3", "8", "2", "5", "4", "7"] ; // vertices balanced, opposed faces imbalance minimised

// for the Octahedron, balancing halves is the same thing as balancing vertices

vertices = [
[0.0, 0.0, C0], //0
[0.0, 0.0,-C0], //1
[ C0,  C1, C2], //2
[ C0, -C1,-C2], //3
[-C1,  C0, C2], //4
[ C1,  C0,-C2], //5
[-C0, -C1, C2], //6
[-C0,  C1,-C2], //7
[ C1, -C0, C2], //8
[-C1, -C0,-C2]  //9
];

faces = [
[ 0 , 8 , 3 , 2 ],
[ 0 , 2 , 5 , 4],
[ 0 , 4 , 7 , 6],
[ 0 , 6 , 9 , 8],
[ 1 , 7 , 4 , 5],
[ 1 , 5 , 2 , 3],
[ 1 , 3 , 8 , 9],
[ 1 , 9 , 6 , 7]
];



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
