// A hull no closed form covers, rendered finely enough that the exporter
// builds it from the children tessellated at its own fineness instead of
// from OpenSCAD's render; the fallback count and the volume are what the
// expected file holds to. Its face count is not dependable to the unit,
// which is why step metrics are compared with slack on that one field.
module trio() {
  sphere(r = 1);
  translate([4, 0, 1]) sphere(r = 1.5);
  translate([1, 3.5, 2]) sphere(r = 0.75);
}
hull($fn = 180) trio($fn = 180);
