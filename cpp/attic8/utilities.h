namespace dnt
{
/// Compute relative orientation of point x wrt segment [a, b]
  double orient1d(double a, double b, double x);

  /// Compute relative orientation of points a, b, c. The orientation
  /// is such that orient2d(a, b, c) > 0 if a, b, c are ordered
  /// counter-clockwise.
  double _orient2d(const double* a, const double* b, const double* c);

  /// Convenience function using dolfin::Point
  double orient2d(const Point& a, const Point& b, const Point& c);

  /// Compute relative orientation of points a, b, c, d. The
  /// orientation is such that orient3d(a, b, c, d) > 0 if a, b, c, d
  /// are oriented according to the left hand rule.
  double _orient3d(const double* a, const double* b, const double* c, const double* d);

  /// Convenience function using dolfin::Point
  double orient3d(const Point& a, const Point& b, const Point& c, const Point& d);
}
