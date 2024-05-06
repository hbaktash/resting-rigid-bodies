// Gyrate Octogonal Dipyramid
// 3 sided faces = 16
// degree 4 vertices = 6
// degree 6 vertices = 4
// halves = 6 (the two gyrated halves + 4 halves based on the the degree 6 vertices)
 
txt_depth = .28;
txt_size = .35;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 27; // scale factor which sets the diameter (distance from one [deg 4]-vertex to the vertex opposite) 
minko = 0.15; //chamfer the edges [0 = disabled]
roll = 0.15; // another effet to round the dice (intersects with a smaller sphere) [disabled if <= 0]
minkfn = 80; //sets the $fn variable for chamfer and the sphere

C2 = 1; // the height of each pyramids in the dipyramid; different from the original coordinates (2)

ordi = C2+minko; //original diameter
scafa = diam*.5/ordi; //scale factor

zint = 0;// z intersect of the ray perpendicular to the faces

//labels = ["1","2","3","4","5","6.","7","8","9.","10","11","12","13","14","15","16"]; //default labeling

// a balanced labelling (halves and vertices)
labels = ["3","12","13","6.","9.","16","7","2","14","5","4","11","8","1","10","15"];


// a labelling where the numbers appear from 1 to 16 in order
//labels = ["1","2","3","4","13","14","15","16","12","11","10","9.","8","7","6.","5"]; //default labeling





vertices = [
[ 0.0, 0.0,  C2],
[ 0.0, 0.0,  -C2],
[ cos(360*0/8),  sin(360*0/8),  0.0], //
[ 0.0 ,  sin(360*1/8), cos(360*1/8) ],
//[ cos(360*1/8),  sin(360*1/8), 0.0 ],
[ cos(360*2/8),  sin(360*2/8),  0.0], //
[ 0.0 ,  sin(360*3/8),  cos(360*3/8)],
[ cos(360*4/8),  sin(360*4/8),  0.0], //
[ cos(360*5/8),  sin(360*5/8),  0.0],
[ cos(360*6/8),  sin(360*6/8),  0.0], 
[ cos(360*7/8),  sin(360*7/8),  0.0]
];
faces = [
[ 0 , 3, 2],
[ 3 , 4, 2],
[ 2 , 4, 5],
[ 2 , 5, 1],
[ 0 , 7, 6],
[ 0 , 8, 7],
[ 0 , 9, 8],
[ 0 , 2, 9],
[ 3 , 0, 6],
[ 4 , 3, 6],
[ 5 , 4, 6],
[ 5 , 6, 1],
[ 6 , 7, 1],
[ 7 , 8, 1],
[ 8 , 9, 1],
[ 9 , 2, 1],
];

//color("orange",0.95) polyhedron(points = vertices, faces=faces,  convexity = 20);

//for(i = [0 : len(vertices) - 1]) translate(vertices[i]) cube(.21);

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
    rotate([180,0,0]) translate(barpof) rotate([0,b,c]) linear_extrude(txt_depth,center=true) text(text=txt,size = txt_size, font=txt_font, halign = "center", valign = "center");
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

