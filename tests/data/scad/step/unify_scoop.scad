// A sphere scooped from a cylinder's end, then halved by a plane through
// the scoop's equator. Unifying the scooped shape once inflated a vertex
// tolerance to half the model, and the plane cut then kept both halves.
difference() {
  difference() {
    rotate(a = 180, v = [1, 0, 1]) cylinder(h = 125, r = 25, center = true);
    translate([-72.5, 0, 0]) sphere(25);
  }
  translate([0, 0, -100]) cube(200, center = true);
}
