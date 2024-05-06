use <dice.scad>

translate([0,0,1])
cut([30,0,-1.5],[180,0,0]) {
    d30();
    union() {
        translate([0,0,-30])
        cube([60,60,60],center=true);
        cuboid([20,20,6],1,center=true);
    }
}

   