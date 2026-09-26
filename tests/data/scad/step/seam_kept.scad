// A tube fused with an oblate ellipsoid, then cut by a ring whose inner
// radius equals the ellipsoid's, so the ellipsoid is trimmed exactly at
// its equator. Splitting that face at its seam produced a STEP file that
// read back 55% too large; the unsplit file reads back exactly.
module ring(id, od, h) { difference() { cylinder(r = od / 2, h = h); translate([0, 0, -0.01]) cylinder(r = id / 2, h = h + 0.02); } }
rotate([0, 180, 0]) difference() {
  union() {
    translate([0, 0, 7.375]) {
      ring(id = 11, od = 35, h = 42.625);
      difference() { scale([1, 1, 0.5]) sphere(r = 14.75); cylinder(r = 5.5, h = 16.75, center = true); }
    }
  }
  ring(id = 29.5, od = 36, h = 48);
}
