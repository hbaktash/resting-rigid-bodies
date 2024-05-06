// Heptagonal Trapezohedron
// 4 sided faces = 18
// degree 3 vertices = 18
// degree 9 vertices = 2

txt_depth = .09;
txt_size = .45;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 30; // scale factor which sets the diameter (distance from one vertex to the vertex opposite) 
minko = 0.1; //rounds the edges
roll = -0.08; // another effect to round the dice (intersects with a smaller sphere) [0 = disabled]
minkfn = 80; //sets the $fn variable for the "rounding"s


C0 = 0.031091204125763378917797589082543290835776502803909279129259;
C1 = 1; // this parameter can be varied to change the shape
C2 = -1.4323768688387176475837417002879003086689156586745832750427;


ordi = C1+minko;
scafa = diam*.5/ordi;

//labels = ["1","2","3","4","5","6.","7","8","9.","10","11","12","13","14","15","16","17","18"]; // default labeling,

// this labeling has balanced opposed faces. Four [degree 3] vertices deviate by 4.5, four deviate by 3.5, two by 2.5 and two by 1.5; the remaing six degree 3 vertices deviate by 0.5, as do the degree 9 vertices (which are also the halves).
labels= ["6.", "3", "9.", "15", "14", "2", "7", "11", "18", "17", "12", "8", "1", "13", "16", "10", "4", "5"];



zint = -(C0+1)/4*C1 - 1/(C2*C1); // z intersect of the perpendicular to a face

vertices = [
[0,0,C1],
[0,0,-C1],
[sin(180/9*1),cos(180/9*1),C0],
[sin(180/9*2),cos(180/9*2),-C0],
[sin(180/9*3),cos(180/9*3),C0],
[sin(180/9*4),cos(180/9*4),-C0],
[sin(180/9*5),cos(180/9*5),C0],
[sin(180/9*6),cos(180/9*6),-C0],
[sin(180/9*7),cos(180/9*7),C0],
[sin(180/9*8),cos(180/9*8),-C0],
[sin(180/9*9),cos(180/9*9),C0],
[sin(180/9*10),cos(180/9*10),-C0],
[sin(180/9*11),cos(180/9*11),C0],
[sin(180/9*12),cos(180/9*12),-C0],
[sin(180/9*13),cos(180/9*13),C0],
[sin(180/9*14),cos(180/9*14),-C0],
[sin(180/9*15),cos(180/9*15),C0],
[sin(180/9*16),cos(180/9*16),-C0],
[sin(180/9*17),cos(180/9*17),C0],
[sin(180/9*18),cos(180/9*18),-C0],
];

faces = [
[ 0 ,  2 ,  3 ,  4 ],
[ 0 ,  4 ,  5 ,  6 ],
[ 0 ,  6 ,  7 ,  8 ],
[ 0 ,  8 ,  9 , 10 ],
[ 0 ,  10 ,  11 , 12 ],
[ 0 ,  12 ,  13 , 14 ],
[ 0 ,  14 ,  15 , 16 ],
[ 0 ,  16 , 17 ,  18 ],
[ 0 ,  18 ,  19 , 2 ],
[ 1 ,  5 ,  4 ,  3],
[ 1 ,  7 ,  6 ,  5],
[ 1 ,  9 ,  8 ,  7],
[ 1 ,  11 , 10 , 9],
[ 1 ,  13 , 12 , 11],
[ 1 ,  15 , 14 , 13],
[ 1 ,  17 , 16 , 15],
[ 1 ,  19 , 18 , 17],
[ 1 , 3 ,  2 , 19]
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