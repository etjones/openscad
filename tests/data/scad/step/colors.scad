color("red") cube(10);
translate([20, 0, 0]) color("blue") sphere(5);
translate([40, 0, 0]) cube(6);
// overlapping: later blue claims the shared material
translate([0, 30, 0]) union() { color("red") cube(10); color("blue") translate([5, 5, 5]) cube(10); }
// preview rule: outer color wins
translate([30, 30, 0]) color("green") color("yellow") cylinder(h = 5, r = 4);
// difference keeps the minuend's color
translate([60, 0, 0]) difference() { color("red") cube(10); color("blue") translate([5, 5, -1]) cylinder(h = 12, r = 2); }
