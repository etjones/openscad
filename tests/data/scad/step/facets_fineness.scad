// Faceting requested through $fa/$fs instead of $fn.
cylinder(d = 9.8, h = 13, $fa = 70);                          // 6 segments
translate([20, 0, 0]) cylinder(r = 3, h = 5, $fs = 3);         // 7 segments
translate([40, 0, 0]) sphere(r = 4, $fa = 45);                 // 8 segments, mesh fallback
translate([60, 0, 0]) cylinder(r = 1, h = 5);                  // default fineness: exact
translate([0, 20, 0]) linear_extrude(3) circle(r = 5, $fa = 60);  // hexagon prism
translate([20, 20, 0]) cylinder(r = 5, h = 5, $fa = 1);        // 16 segments: exact
