// Biaugmented Hexagonal Antiprism
// 3 sided faces = 24
// degree 5 vertices = 12
// degree 6 vertices = 2

txt_depth = .15;
txt_size = .4;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 24; // scale factor which sets the diameter (distance from one vertex to the vertex opposite) 
minko = 0.1; //rounds the edges
roll = -0.15; // another effet to round the dice (intersects with a smaller sphere) [0 = disabled]
minkfn = 80; //sets the $fn variable for the "rounding"s



C0 = 1.6528916502810694800982406465853878177773;
C1 = 1.0478913165750138709807423198640911228892;

ordi = 2*C0+minko;
scafa = diam*.5/ordi;

labels = ["1","2","3","4","5","6.","7","8","9.","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24"]; // default labeling

zint = 0; //-(C1 - (1+2*cos(36))/C1)/4; // z intersect of the perpendicular to a face

vertices = [
[0,0,C0-C1/2],
[0,0,C1/2-C0],
[sin(30),cos(30),C1/2],
[sin(60),cos(60),-C1/2],
[sin(90),cos(90),C1/2],
[sin(120),cos(120),-C1/2],
[sin(150),cos(150),C1/2],
[sin(180),cos(180),-C1/2],
[sin(210),cos(210),C1/2],
[sin(240),cos(240),-C1/2],
[sin(270),cos(270),C1/2],
[sin(300),cos(300),-C1/2],
[sin(330),cos(330),C1/2],
[sin(360),cos(360),-C1/2],
];

faces = [
[0 , 2 , 4],
[0 , 4 , 6],
[0 , 6 , 8],
[0 , 8 , 10],
[0 , 10 , 12],
[0 , 12 , 2],
[1 , 5 , 3],
[1 , 7 , 5],
[1 , 9 , 7],
[1 , 11 , 9],
[1 , 13 , 11],
[1 , 3 , 13],
[2 , 3 , 4],
[5 , 4 , 3],
[4 , 5 , 6],
[7 , 6 , 5],
[6 , 7 , 8],
[9 , 8 , 7],
[8 , 9 , 10],
[11 , 10 , 9],
[10 , 11 , 12],
[13 , 12 , 11],
[12 , 13 , 2],
[3 , 2 , 13]
];

//color("orange",0.5) 
//polyhedron(points = vertices, faces=faces,  convexity = 20);

function add3(v, i = 0, r) = 
    i < len(v) ? 
        i == 0 ?
            add3(v, 1, v[0]) :
            add3(v, i + 1, r + v[i]) :
        r;

function facecoord(n) = [for(i = [0 : len(faces[n]) - 1]) vertices[faces[n][i]] ];

//echo(add3(facecoord(4))/4);

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
    for(i=[0:len(faces)-1]) facetext(facecoord(i),labels[i],faces[i][0]);
}