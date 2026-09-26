// Three touching bodies: a revolved ring, a plate on its outside, and a
// wedge sharing the plate's face. OCCT's plain fuse silently dropped the
// plate here (its volume stayed within the bounds the inputs imply and
// the piece count did not rise), so the union is also checked for
// operands that are no longer inside the result.
ringRadius = 5.72171;
union() {
  translate([0, 0, 9.5]) rotate_extrude() translate([ringRadius, 0]) circle(r = 0.5);
  translate([ringRadius, -1, 4.5]) cube([0.5, 2, 5.5]);
  translate([ringRadius, -0.5, 4.5]) rotate([90, 0, 180]) linear_extrude(1)
    polygon([[0, 0], [0, 5.5], [1, 0]]);
}
