hull() for (x = [0, 30], y = [0, 20]) translate([x, y, 0]) { cylinder(h = 8, r = 3, $fn = 64); translate([0, 0, 8]) sphere(3, $fn = 64); }
