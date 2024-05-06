phi = (1+sqrt(5))/2;

module label1(v,r,sz) {
  rotate(v[0])
    translate([r-1,0,0])
    rotate([0,90,0])
    linear_extrude(1)
    text(v[1],valign="center",halign="center",font="Times New Roman:style=Bold",size=sz);
}

module label2(v,r,sz) {
  rotate(v[0])
    for(j=[0:len(v[1])-1]) {
      translate([r-1,0,0])
	rotate([-j*360/len(v[1]),0,0])
	rotate([0,90,0])
	translate([0,10,0])
	linear_extrude(1)
	text(v[1][j],valign="bottom",halign="center",font="Times New Roman:style=Bold",size=sz);
    }
}

function auto_labels(n) =
  [for(i=[1:n]) (n<9 || (i!=6 && i!=9))?str(i):str(i,".")];

module gen_die(vv,height=50,mn=0,tsize=0,s=35,simp=true) {
  r = height/2;
  m = (mn==0)?height/20:mn;
  sz = (tsize==0)?height/5:tsize;
  difference() {
    minkowski() {
      intersection() {
	intersection_for(v=vv) {
	  rotate(v[0])
	    translate([-2*r-m,0,0])
	    cube([6*r,6*r,6*r],center=true);
	}
	sphere(s,$fa=1);
      }
      sphere(m,$fn=30);
    }
    for(v=vv) {
      if (simp) {
	label1(v,r,sz);
      } else {
	label2(v,r,sz);
      }
    }
  }
}

module cuboid(dim, r, center=false) {
  minkowski() {
    d=[dim[0]-2*r,dim[1]-2*r,dim[2]-2*r];
    cube(d, center=center);
    sphere(r, $fn=30);
  }
}

module cut(t,r) {
  translate(t)
    difference() {
    children(0);
    children(1);
  }
  translate(-t)
    rotate(r)
    intersection() {
    children(0);
    children(1);
  }
}

function zip(a, b) =
  [for (i=[0:len(a)-1]) [a[i], b[i]]];

d4_labels = [["1","2","3"],["4","3","2"],["4","2","1"],["4","1","3"]];
    
module d4(size = 50, labels=d4_labels) {
  a=asin(sqrt(1/3));
  s=size/78;
  angles = [
	    [88,a,45],[88,a,225],
	    [32,-a,135],[152,-a,315]
	    ];
  v = zip(angles,labels);

  scale([s,s,s])
    gen_die(v,tsize=20,simp=false,s=50);
}

d4b_labels = auto_labels(4);

module d4b(size=50, labels=d4b_labels, tsize=7) {
  a=asin(sqrt(1/3));
  height=size;
    
  angles = [
	    [88,a,45],[88,a,225],
	    [32,-a,135],[152,-a,315]
	    ];
  vv = zip(angles,labels);
  r = 9*height/20;
  m = height/20;
  sz = tsize;
  difference() {
    minkowski() {
      intersection_for(v=vv) {
	rotate(v[0])
	  translate([r/3,0,0])
	  cube([2*r,6*r,6*r],center=true);
      }
      sphere(m,$fn=30);
    }
    for(v=vv) {
      label1(v,4*r/3+m,sz);
    }
  }    
}

d6_labels = auto_labels(6);

module d6(size=50, labels=d6_labels, tsize=20) {
  s=size/55;

  angles = [[0,0,0],[90,0,90],[0,90,90],
	    [180,270,90],[90,180,90],[0,0,180] ];
    
  desc = zip(angles, labels);
  rotate([45,asin(sqrt(1/3)),0])
    scale([s,s,s])
    gen_die(desc,height=50,tsize=tsize);
}

d8_labels = auto_labels(8);

module d8(size=50, labels=d8_labels, tsize=15) {
  s1=asin(sqrt(1/3));
  s2=-s1;

  angles = [[30,s1,45],[30,s1,135],[30,s1,225],[30,s1,315],
	    [-30,s2,135],[-30,s2,45],[-30,s2,315],[-30,s2,225]];

  desc = zip(angles, labels);
  s = size/55;
  scale([s,s,s])
    gen_die(desc,tsize=tsize,height=50);
}

function bipyramidal(n, th=40) =
  [for(i=[0:n-1])
      ((i%2)==0)?
	[0,th,360*i/n]:
	[0,-th,360*(i-1)/n]];
 
function trapezohedron(n, th=40, rot=90) =
  [for(i=[1:n])
      (2*i<=n)?
	[-rot,th,720*i/n]:
	[rot,-th,360*(2*(n-i)-3)/n]];

module trap_die(n, th=40, rot=0, size=50, labels=[], tsize=13) {
  angles = trapezohedron(n,th,rot);
  l = (len(labels)==0)?auto_labels(n):labels;
  desc = zip(angles, l);
  s=size/55;
  scale([s,s,s])
    gen_die(desc,tsize=tsize,s=32);
}

d10_labels = auto_labels(10);

module d10(size=50, labels=d10_labels, tsize=13) {
  trap_die(10, size=size, labels=labels, tsize=tsize);
}

d12_labels = auto_labels(12);

module d12(size=50, labels=d12_labels, tsize=12) {
  s1=-atan(1/2);
  s2=-s1;

  angles = [[0,-90,0],[0,s1,0],[0,s1,72],[0,s1,144],
	    [0,s1,216],[0,s1,288],[0,s2,108],[0,s2,36],
	    [0,s2,324],[0,s2,252],[0,s2,180],[0,90,0]];

  desc = zip(angles, labels);
  s=size/55;
  scale([s,s,s])
    rotate([0,38,0])
    gen_die(desc,tsize=tsize,s=27);
}

d20_labels = auto_labels(20);

module d20(size=50, labels=d20_labels, tsize=9) {
  l1=atan((3+sqrt(5))/4);
  l2=atan((3-sqrt(5))/4);

  angles = [[90,-l1,0],[90,-l1,72],[90,-l1,144],[90,-l1,216],
	    [90,-l1,288],[-90,-l2,0],[-90,-l2,72],[-90,-l2,144],
	    [-90,-l2,216],[-90,-l2,288],[90,l2,108],[90,l2,36],
	    [90,l2,324],[90,l2,252],[90,l2,180],[-90,l1,108],
	    [-90,l1,36],[-90,l1,324],[-90,l1,252],[-90,l1,180]];

  desc = zip(angles, labels);    
  s=size/55;
  scale([s,s,s])
    gen_die(desc,tsize=tsize,s=27);
}

d24_labels = auto_labels(24);

module d24(size=50,labels=d24_labels) {
  l1=atan(1/3);
  l2=90-atan(1/2);
  angles=[[0,l2,-135],[0,l2,-45],[0,l2,45],[0,l2,135],
	  [0,l1,-90-l1],[0,l1,-90+l1],[0,l1,-l1],[0,l1,l1],
	  [0,l1,90-l1],[0,l1,90+l1],[0,l1,180-l1],[0,l1,180+l1],
	  [0,-l1,l1],[0,-l1,-l1],[0,-l1,-90+l1],[0,-l1,-90-l1],
	  [0,-l1,180+l1],[0,-l1,180-l1],[0,-l1,90+l1],[0,-l1,90-l1],
	  [0,-l2,-45],[0,-l2,-135],[0,-l2,135],[0,-l2,45]];
  desc = zip(angles,labels);
  s=size/55;
  rotate([45,0,0])
    scale([s,s,s])
    gen_die(desc,tsize=7,s=25);
}

d30_labels = auto_labels(30);

module d30(size=50, labels=d30_labels) {
  s1=atan(phi);
  s2=atan(1/phi);
  s3=0;
  s4=180;
  s5=-s2;
  s6=-s1;

  angles = [[0,s1,0],[0,s1,72],[0,s1,144],[0,s1,216],
	    [0,s1,288],[90,s2,36],[90,s2,108],[90,s2,180],
	    [90,s2,252],[90,s2,324],[s1+90,s3,18],[s1+90,s3,90],
	    [s1+90,s3,162],[s1+90,s3,234],[s1+90,s3,306],[s2,s4,306],
	    [s2,s4,234],[s2,s4,162],[s2,s4,90],[s2,s4,18],
	    [90,s5,144],[90,s5,72],[90,s5,0],[90,s5,288],
	    [90,s5,216],[0,s6,108],[0,s6,36],[0,s6,324],
	    [0,s6,252],[0,s6,180]];

  desc = zip(angles, labels);

  s=size/55;
  scale([s,s,s])
    gen_die(desc,tsize=7,s=25);
}

module support(z) {
  translate([0,0,-z-5]) {
    cylinder(1,25,25);
    cylinder(z+5,1,1,$fn=12);
    cylinder(5,25,2);
  }
}

function rot(v) =
  [[v[0],v[1],v[2]],
   [v[1],v[2],v[0]],
   [v[2],v[0],v[1]]];

function c2a(v) =
  [0,
   atan2(v[2],sqrt(v[0]*v[0]+v[1]*v[1])),
   atan2(v[1],v[0])];

function strip(v) = len(v)==1?v[0]:
  concat(strip([for(i=[0:len(v)-2]) v[i]]),v[len(v)-1]);
       
function signs(v) =
  [[v[0],v[1],v[2]],
   [v[0],v[1],-v[2]],
   [v[0],-v[1],v[2]],
   [v[0],-v[1],-v[2]],
   [-v[0],v[1],v[2]],
   [-v[0],v[1],-v[2]],
   [-v[0],-v[1],v[2]],
   [-v[0],-v[1],-v[2]]];

d48_labels = auto_labels(48);

module d48(size=50, labels=d48_labels) {
  v1=strip([for (i=rot([1,1+sqrt(2),1+2*sqrt(2)])) signs(i)]);
  v2=strip([for (i=rot([1,1+2*sqrt(2),1+sqrt(2)])) signs(i)]);

  v=strip([v1,v2]);

  desc = zip([for(i=v) c2a(i)], labels);
    
  s=size/50;
  scale([s,s,s])
    gen_die(desc,tsize=4, s=25);
}

d60_labels = auto_labels(60);

module d60(size=50, labels=d60_labels) {
  v0=[[0,1,3*phi],[0,1,-3*phi],
      [0,-1,3*phi],[0,-1,-3*phi]];
  v1=strip([for(i=v0) rot(i)]);
  v2=strip([for (i=rot([1,2+phi,2*phi])) signs(i)]);
  v3=strip([for (i=rot([phi,2,2*phi+1])) signs(i)]);

  v=strip([v1,v2,v3]);

  desc = zip([for(i=v) c2a(i)], labels);
    
  s=size/50;
  scale([s,s,s])
    gen_die(desc,tsize=4,s=24);
}

d120_labels = auto_labels(120);

module d120(size=50, labels=d120_labels) {
  v1=strip([for (i=rot([1/phi,1/phi,3+phi])) signs(i)]);
  v2=strip([for (i=rot([2/phi,phi,1+2*phi])) signs(i)]);
  v3=strip([for (i=rot([1/phi,phi*phi,-1+3*phi])) signs(i)]);
  v4=strip([for (i=rot([-1+2*phi,2,2+phi])) signs(i)]);
  v5=strip([for (i=rot([phi,3,2*phi])) signs(i)]);

  v=strip([v1,v2,v3,v4,v5]);

  desc = zip([for(i=v) c2a(i)], labels);
    
  s=size/50;
  scale([s,s,s])
    gen_die(desc,tsize=2.5,s=24);
}

module all() {
  dist = 100;
  translate([dist,0,0]) d4();
  translate([dist/2,dist*sqrt(3/4),0]) d4b();
  translate([-dist/2,dist*sqrt(3/4),0]) d6();
  translate([-dist,0,0]) d8();
  translate([-dist/2,-dist*sqrt(3/4),0]) d10();
  translate([dist/2,-dist*sqrt(3/4),0]) d12();
  translate([dist/2,dist*sqrt(1/12),dist*sqrt(2/3)]) d20();
  translate([-dist/2,dist*sqrt(1/12),dist*sqrt(2/3)]) d24();
  translate([0,-dist*sqrt(1/3),dist*sqrt(2/3)]) d30();
  translate([0,dist*sqrt(1/3),-dist*sqrt(2/3)]) d48();
  translate([-dist/2,-dist*sqrt(1/12),-dist*sqrt(2/3)]) d60();
  translate([dist/2,-dist*sqrt(1/12),-dist*sqrt(2/3)]) d120();
}

d100_labels = auto_labels(100);

module d100(size=50, labels=d100_labels) {
    v=[1,6,11,15,17,17,15,11,6,1];
    function ang(n,a) = [for(i=[0:n-1])[0,a,i*360/n]];
    angles=strip([for(i=[0:len(v)-1]) ang(v[i],90-360*i/18)]);
    desc = zip(angles,labels);
    s=size/50;
    scale([s,s,s])
        gen_die(desc,tsize=4,s=23);
};