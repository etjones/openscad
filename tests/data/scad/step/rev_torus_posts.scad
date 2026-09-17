hull() for (x = [0, 30], y = [0, 20]) translate([x, y, 0]) minkowski() { cylinder(h = 6, r = 2, $fn = 64); sphere(2, $fn = 64); }
