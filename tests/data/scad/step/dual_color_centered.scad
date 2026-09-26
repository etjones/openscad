color("red") {
  difference() {
    cube(10);
    cylinder(h=15, r=5);
  }
}

color("blue") {
  translate([0, 0, 0]) cylinder(h=15, r=5);
}

translate([20, 0, 0]) {
  union() {
    color("red") {
      cube(10, center=true);
    }

    color("magenta") {
      sphere(6);
    }
  }
}
  