// Octogonal Dipyramid
// 3 sided faces = 16
// degree 4 vertices = 8
// degree 8 vertices = 2
// halves = 10 (two are the degree 8 vertices)
 
txt_depth = .11;
txt_size = .35;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 30; // scale factor which sets the diameter (distance from one [deg 4]-vertex to the vertex opposite) 
minko = 0.15; //chamfer the edges [0 = disabled]
roll = 0.35; // another effet to round the dice (intersects with a smaller sphere) [disabled if <= 0]
minkfn = 80; //sets the $fn variable for chamfer and the sphere

C0 = 1/sqrt(3);
C1 = 2/sqrt(3);
C2 = 1; // the height of each pyramids in the dipyramid; different from the original coordinates (2)

ordi = max(C2,C1)+minko; //original diameter
scafa = diam*.5/ordi; //scale factor

zint = 0; //-C2/3 + 2/(3*C2);//1/sqrt(27); // z intersect of the ray perpendicular to the faces

//labels = ["1","2","3","4","5","6.","7","8","9.","10","11","12","13","14","15","16"]; //default labeling

//This labelling has balanced opposed faces. Two vertices devaite by 2 points, and four halves deviate by 2 points.
labels = ["13","12","8","1","14","11","7","2","3","6.","10","15","4","5","9.","16"];


//This labelling has balanced vertices and halves. The opposed faces are not balanced (they all deviate by 1 point).
//labels = ["7","1","11","14","8","2","12","13","10","16","6.","3","9.","15","5","4"];


vertices = [
[ 0.0, 0.0,  C2],
[ 0.0, 0.0,  -C2],
[ cos(360*1/8),  sin(360*1/8),  0.0],
[ cos(360*2/8),  sin(360*2/8),  0.0],
[ cos(360*3/8),  sin(360*3/8),  0.0],
[ cos(360*4/8),  sin(360*4/8),  0.0],
[ cos(360*5/8),  sin(360*5/8),  0.0],
[ cos(360*6/8),  sin(360*6/8),  0.0],
[ cos(360*7/8),  sin(360*7/8),  0.0],
[ cos(360*8/8),  sin(360*8/8),  0.0]
];
faces = [
[ 0 , 3, 2],
[ 0 , 4, 3],
[ 0 , 5, 4],
[ 0 , 6, 5],
[ 0 , 7, 6],
[ 0 , 8, 7],
[ 0 , 9, 8],
[ 0 , 2, 9],
[ 2 , 3, 1],
[ 3 , 4, 1],
[ 4 , 5, 1],
[ 5 , 6, 1],
[ 6 , 7, 1],
[ 7 , 8, 1],
[ 8 , 9, 1],
[ 9 , 2, 1],
];



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

