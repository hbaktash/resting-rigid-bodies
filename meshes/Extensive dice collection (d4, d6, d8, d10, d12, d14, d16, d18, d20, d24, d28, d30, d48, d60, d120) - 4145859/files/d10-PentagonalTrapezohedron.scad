// Pentagonal Trapezohedron
// 4 sided faces = 10
// degree 3 vertices = 10
// degree 5 vertices = 2

txt_depth = .09;
txt_size = .4;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 24; // scale factor which sets the diameter (distance from one vertex to the vertex opposite) 
minko = 0.1; //rounds the edges
roll = 0.15; // another effet to round the dice (intersects with a smaller sphere) [0 = disabled]
minkfn = 80; //sets the $fn variable for the "rounding"s


C0 = 0.1067764671821527129104854331344823715423;

ordi = 1+C0*C0+minko;
scafa = diam*.5/ordi;

//labels = ["1","2","3","4","5","6.","7","8","9.","10"]; // default labeling

//the upcoming labeling has balanced opposed faces and the smallest imbalance of the vertices (.5 for the degree 5 vertices; two degree 3 vertices deviate by 2.5, two by 1.5 and the remaining six by .5)

labels = ["1","8","5","4","9.","7","2","10","3","6."];

C1 = 1+C0+C0*C0;
zint = -(C1 - (1+2*cos(36))/C1)/4; // z intersect of the perpendicular to a face

vertices = [
[0,0,1+C0*C0],
[0,0,-(1+C0*C0)],
[sin(36),cos(36),C0],
[sin(72),cos(72),-C0],
[sin(108),cos(108),C0],
[sin(144),cos(144),-C0],
[sin(180),cos(180),C0],
[sin(216),cos(216),-C0],
[sin(252),cos(252),C0],
[sin(288),cos(288),-C0],
[sin(324),cos(324),C0],
[sin(360),cos(360),-C0],
];

faces = [
[ 0 ,  2 ,  3 ,  4 ],
[ 0 ,  4 ,  5 ,  6 ],
[ 0 ,  6 ,  7 ,  8 ],
[ 0 ,  8 ,  9 , 10 ],
[ 0 , 10 , 11 ,  2 ],
[ 1 ,  5 ,  4 ,  3],
[ 1 ,  7 ,  6 ,  5],
[ 1 ,  9 ,  8 ,  7],
[ 1 ,  11 , 10 , 9],
[ 1 , 3 ,  2 , 11]
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