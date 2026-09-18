// A non-uniformly scaled sphere fused to a cylinder. The ellipsoid's face
// wrapped all the way around its own seam, and the STEP file written from
// it read back with the ellipsoid empty. Faces are split at the seam
// before writing, so built and roundtrip below must both hold it.
union() {
  translate([0, 0, -11]) rotate([0, 90, 0]) scale([0.5, 1, 1]) sphere(r = 2);
  translate([0, 0, -15]) cylinder(h = 4, r = 1.2);
}
