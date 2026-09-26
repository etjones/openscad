hull() for (x = [0, 30], y = [0, 30]) translate([x, y, 0]) { cylinder(h = 5, r = 4, $fn = 64); translate([0, 0, 5]) cylinder(h = 3, r1 = 4, r2 = 2, $fn = 64); }
