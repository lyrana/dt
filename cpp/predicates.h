// This is a DnT header file for predicates.cpp which provides
//
//   Routines for Arbitrary Precision Floating-point Arithmetic
//   and Fast Robust Geometric Predicates
//
// by
//
//   Jonathan Richard Shewchuk
//
// Code is placed in the public domain.

#ifndef DNT_PREDICATES_H
#define DNT_PREDICATES_H

namespace dnt
{
  /// Compute relative orientation of point x wrt segment [a, b]
  double orient1d(double a, double b, double x);

  /// Compute relative orientation of points a, b, c. The orientation
  /// is such that orient2d(a, b, c) > 0 if a, b, c are ordered
  /// counter-clockwise.
  double _orient2d(const double* a, const double* b, const double* c);

  /// Compute relative orientation of points a, b, c, d. The
  /// orientation is such that orient3d(a, b, c, d) > 0 if a, b, c, d
  /// are oriented according to the left hand rule.
  double _orient3d(const double* a, const double* b, const double* c, const double* d);

  // Implementation of collision detection predicates
  
  bool _collides_segment_point_1d(double p0,
                                  double p1,
                                  double point);
}

namespace dolfin
{

  /// Initialize tolerances for exact arithmetic
  void exactinit();

  /// Class used for automatic initialization of tolerances at startup.
  /// A global instance is defined inside predicates.cpp to ensure that
  /// the constructor and thus exactinit() is called.

  class PredicateInitialization
  {
  public:

    PredicateInitialization()
    {
      exactinit();
    }

  };

}

#endif
