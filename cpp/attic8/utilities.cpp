// /home/tph/workspace/dolfin/dolfin/geometry/predicates.cpp


//-----------------------------------------------------------------------------
double dnt::orient1d(double a, double b, double x)
{
  if (x > std::max(a, b)) return 1.0;
  if (x < std::min(a, b)) return -1.0;
  return 0.0;
}
//-----------------------------------------------------------------------------
double dnt::orient2d(const Point& a, const Point& b, const Point& c)
{
  return dnt::_orient2d(a.coordinates(),
                           b.coordinates(),
                           c.coordinates());
}
//-----------------------------------------------------------------------------
double dnt::orient3d(const Point& a, const Point& b, const Point& c, const Point& d)
{
  return dnt::_orient3d(a.coordinates(),
                           b.coordinates(),
                           c.coordinates(),
                           d.coordinates());
}


REAL dnt::_orient2d(const REAL *pa, const REAL *pb, const REAL *pc)
/* REAL *pa; */
/* REAL *pb; */
/* REAL *pc; */
{
  REAL detleft, detright, det;
  REAL detsum, errbound;

  detleft = (pa[0] - pc[0]) * (pb[1] - pc[1]);
  detright = (pa[1] - pc[1]) * (pb[0] - pc[0]);
  det = detleft - detright;

  if (detleft > 0.0) {
    if (detright <= 0.0) {
      return det;
    } else {
      detsum = detleft + detright;
    }
  } else if (detleft < 0.0) {
    if (detright >= 0.0) {
      return det;
    } else {
      detsum = -detleft - detright;
    }
  } else {
    return det;
  }

  errbound = ccwerrboundA * detsum;
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  return orient2dadapt(pa, pb, pc, detsum);
}
