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


C0 = 0.07217073141842624582587810305031308204629;

ordi = 1+C0*C0+minko;
scafa = diam*.5/ordi;

//labels = ["1","2","3","4","5","6.","7","8","9.","10","11","12"]; // default labelling

//Here is a labelling with balanced vertices of degree 5; vertices of degree 3 deviate by 1.5 (four times) or .5 (eight times); halves deviate by 2 (four times), 1 (six times) or 0 (twice).
labels = ["1","7","12","4","6.","9.","11","2","5","10","3","8"];

//The following labelling tries to minimises the unbalance of the halves (the way they are defined here, they cannot be balanced): four are balanced, eight deviate by 1 point; vertices of degree 5 deviate by 1 point; vertices if degree 3 deviate by 3.5 (once), 2.5 (thrice), 1.5 (thrice) and .5 (five times).
//labels = ["1","7","9.","3","8","10","12","4","5","11","2","6."];

C1 = 1+C0+C0*C0;
zint = -(C1 - (1+2*cos(36))/C1)/4; // z intersect of the perpendicular to a face

vertices = [
[0,0,1+C0*C0],
[0,0,-(1+C0*C0)],
[sin(30),cos(30),C0],
[sin(60),cos(60),-C0],
[sin(90),cos(90),C0],
[sin(120),cos(120),-C0],
[sin(150),cos(150),C0],
[sin(180),cos(180),-C0],
[sin(210),cos(210),C0],
[sin(240),cos(240),-C0],
[sin(270),cos(270),C0],
[sin(300),cos(300),-C0],
[sin(330),cos(330),C0],
[sin(360),cos(360),-C0],
];

faces = [
[ 0 ,  2 ,  3 ,  4 ],
[ 0 ,  4 ,  5 ,  6 ],
[ 0 ,  6 ,  7 ,  8 ],
[ 0 ,  8 ,  9 , 10 ],
[ 0 , 10 , 11 ,  12 ],
[ 0 , 12 , 13 ,  2 ],
[ 1 ,  5 ,  4 ,  3],
[ 1 ,  7 ,  6 ,  5],
[ 1 ,  9 ,  8 ,  7],
[ 1 ,  11 , 10 , 9],
[ 1 ,  13 , 12 , 11],
[ 1 , 3 ,  2 , 13]
];


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
    barpof = [for(j=[0:2]) barpo[j]*lenfac*-1]; // stretching coordiante to compensate for chamfering
    translate(barpof) rotate([180,b,c]) linear_extrude(txt_depth,center=false) text(text=txt,size = txt_size, font=txt_font, halign = "center", valign = "center");
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