// Triangular Dipyramid
// 3 sided faces = 6
// degree 4 vertices = 3
// degree 3 vertices = 2
// halves = degree 3 vertices
 
txt_depth = .2;
txt_size = .55;
txt_font = "Arial:style=Bold";
// Warning: a scale factor is applied later (see below)
diam = 22; // scale factor which sets the diameter (distance from one [deg 4]-vertex to the vertex opposite) 
minko = 0.5; //chamfer the edges [0 = disabled]
roll = 0.25; // another effet to round the dice (intersects with a smaller sphere) [disabled if <= 0]
minkfn = 80; //sets the $fn variable for chamfer and the sphere

C2 = 1; // the height of each pyramids in the dipyramid

ordi = max(C2,1)+minko; //original diameter
scafa = diam*.5/ordi; //scale factor

zint = -C2/2 + 1/2; // z intersect of the ray perpendicular to the faces

//labels = ["1","2","3","4","5","6"]; //default labeling

//this labelling minimises the imbalance at the degree 3 vertices; degree 4 vertices are balanced
labels = ["1","5","4","6","2","3"];

vertices = [
[ 0.0, 0.0,  C2],
[ 0.0, 0.0,  -C2],
[ cos(0), sin(0),  0.0],
[ cos(120), sin(120),  0.0],
[ cos(240), sin(240),  0.0]];
faces = [
[ 3 , 2, 0],
[ 4 , 3, 0],
[ 2 , 4, 0],
[ 2 , 3, 1],
[ 3 , 4, 1],
[ 4 , 2, 1]];



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
    lenfac = (norm([.5,C2*.5])+minko)/norm(barpo);
    barpof = [for(j=[0:2]) barpo[j]*lenfac*-1]; // stretching coordiante to compensate for chamfering
    translate(barpof) rotate([180,b,c]) linear_extrude(txt_depth,center=true) text(text=txt,size = txt_size, font=txt_font, halign = "center", valign = "center");
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
    for(i=[0:len(faces)-1]) 
        facetext(facecoord(i),labels[i],faces[i][2]);
}
