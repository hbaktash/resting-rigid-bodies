// Biaugmented Square Antiprism
// 3 sided faces = 16
// degree 5 vertices = 8
// degree 4 vertices = 2

txt_depth = .25;
txt_size = .4;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 24; // scale factor which sets the diameter (distance from one vertex to the vertex opposite) 
minko = 0.1; //rounds the edges
roll = -0.15; // another effet to round the dice (intersects with a smaller sphere) [0 = disabled]
minkfn = 80; //sets the $fn variable for the "rounding"s



C0 = 1.5537739740300373073441589530631469481646;
C1 = 0.9101797211244546826087155156449371392404;

ordi = C0+minko;
scafa = diam*.5/ordi;

//labels = ["1","2","3","4","5","6.","7","8","9.","10","11","12","13","14","15","16"]; // default labeling

zint = 0; //-(C1 - (1+2*cos(36))/C1)/4; // z intersect of the perpendicular to a face

// A labeling with balanced opposed faces and vertices (up to the unavoidable deviation by 0.5).
labels = ["16", "7", "1", "10", "13", "2", "4", "15", "6.", "5", "9.", "14", "11", "12", "8", "3"];


vertices = [
[0,0,C0-C1/2],
[0,0,C1/2-C0],
[sin(45),cos(45),C1/2],
[sin(90),cos(90),-C1/2],
[sin(135),cos(135),C1/2],
[sin(180),cos(180),-C1/2],
[sin(225),cos(225),C1/2],
[sin(270),cos(270),-C1/2],
[sin(315),cos(315),C1/2],
[sin(360),cos(360),-C1/2],
];

faces = [
[0 , 2 , 4],
[0 , 4 , 6],
[0 , 6 , 8],
[0 , 8 , 2],
[1 , 5 , 3],
[1 , 7 , 5],
[1 , 9 , 7],
[1 , 3 , 9],
[2 , 3 , 4],
[5 , 4 , 3],
[4 , 5 , 6],
[7 , 6 , 5],
[6 , 7 , 8],
[9 , 8 , 7],
[8 , 9 , 2],
[3 , 2 , 9]
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
    bar = add3(vert)/len(vert); // barycentre
    length = norm(bar);     // radial distance
    b = acos(bar.z/length); // inclination angle
    c = atan2(bar.y,bar.x); // azimuthal angle
    lenfac = (norm(bar)+minko)/norm(bar);
    barpof = [for(j=[0:2]) bar[j]*lenfac]; // stretching coordinate to compensate for chamfering
    rotate([180,0,0]) 
    translate(barpof) rotate([0,b,c]) linear_extrude(txt_depth,center=false) text(text=txt,size = txt_size, font=txt_font, halign = "center", valign = "center");
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