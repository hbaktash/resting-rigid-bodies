// Heptagonal Trapezohedron
// 4 sided faces = 14
// degree 3 vertices = 14
// degree 7 vertices = 2

txt_depth = .09;
txt_size = .45;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 30; // scale factor which sets the diameter (distance from one vertex to the vertex opposite) 
minko = 0.15; //rounds the edges
roll = -0.08; // another effet to round the dice (intersects with a smaller sphere) [0 = disabled]
minkfn = 80; //sets the $fn variable for the "rounding"s


C0 = 0.0520950836016870306476547719793196081537;
C1 = 1; // this parameter can be varied to change the shape
C2 = -1.5019535518686029797197730724235309017711;


ordi = C1+minko;
scafa = diam*.5/ordi;

//labels = ["1","2","3","4","5","6.","7","8","9.","10","11","12","13","14"]; // default labeling

// this labeling has balanced opposed faces. Four [degree 3] vertices deviate by 2.5 and four deviate by 1.5; the rest, including the degree 7 vertices, deviate by 0.5.
labels = ["2", "10", "11", "3", "6.", "14", "7", "9.", "1", "8", "13", "5", "4", "12"];

zint = -(C0+1)/4*C1 - 1/(C2*C1); // z intersect of the perpendicular to a face

vertices = [
[0,0,C1],
[0,0,-C1],
[sin(180/7*1),cos(180/7*1),C0],
[sin(180/7*2),cos(180/7*2),-C0],
[sin(180/7*3),cos(180/7*3),C0],
[sin(180/7*4),cos(180/7*4),-C0],
[sin(180/7*5),cos(180/7*5),C0],
[sin(180/7*6),cos(180/7*6),-C0],
[sin(180/7*7),cos(180/7*7),C0],
[sin(180/7*8),cos(180/7*8),-C0],
[sin(180/7*9),cos(180/7*9),C0],
[sin(180/7*10),cos(180/7*10),-C0],
[sin(180/7*11),cos(180/7*11),C0],
[sin(180/7*12),cos(180/7*12),-C0],
[sin(180/7*13),cos(180/7*13),C0],
[sin(180/7*14),cos(180/7*14),-C0]
];

faces = [
[ 0 ,  2 ,  3 ,  4 ],
[ 0 ,  4 ,  5 ,  6 ],
[ 0 ,  6 ,  7 ,  8 ],
[ 0 ,  8 ,  9 , 10 ],
[ 0 ,  10 ,  11 , 12 ],
[ 0 ,  12 ,  13 , 14 ],
[ 0 , 14 , 15 ,  2 ],
[ 1 ,  5 ,  4 ,  3],
[ 1 ,  7 ,  6 ,  5],
[ 1 ,  9 ,  8 ,  7],
[ 1 ,  11 , 10 , 9],
[ 1 ,  13 , 12 , 11],
[ 1 ,  15 , 14 , 13],
[ 1 , 3 ,  2 , 15]
];

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