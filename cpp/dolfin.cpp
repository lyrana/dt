/*! \file dolfin.cpp

  \brief This file has the source code for a C++ implementation of dolfin-based functions.

  The following C++ functions are defined:
    divide_by_cell_volumes()
    interpolate_field_to_points()
    cell_contains_point()
    particle_stays_in_cell()

  \namespace dnt
  \sa dolfin.h, dolfin_cpp.cpp, 

*/

#include "dolfin.h"
#include "predicates.h"

namespace py = pybind11;

/*
void init_ex2(py::module &m) {
m.def("sub", [](int a, int b) { return a - b; });
}
*/

namespace dnt
{

  //! Test whether a point lies within the cell defined by the vertices.
  /*!
    Overload on the number of particle coordinates: 1D, 2D, 3D 
  */

  // v1:           cell_contains_point(mesh, vertices,  psegOut[ipOut]);


  bool cell_contains_point(void)
  {
    return true;
  }
  
  bool cell_contains_point(dolfin::Mesh& mesh,
                           const unsigned int* vertices,
                           double* point)
  {
    // The algorithm used depends on the mesh dimension
    const std::size_t tdim = mesh.topology().dim();
    const std::size_t gdim = mesh.geometry().dim();

    //v1:    auto v = vertices.unchecked<1>(); // <1> is the number of array indices. See pybind11 docs.

    // 1D cell
    if (tdim == 1 && gdim == 1)
      {
        //v1:        auto v0 = v(0);
        //v1:        auto v1 = v(1);
        auto v0 = vertices[0];
        auto v1 = vertices[1];
        double p0 = mesh.coordinates()[v0];
        double p1 = mesh.coordinates()[v1];

        std::cout << "point[0] = " << point[0] << " p0 = " << p0 << ", p1 = " << p1 << std::endl;
    
        
        if (p0 > p1)
          std::swap(p0, p1);
        return p0 <= point[0] and point[0] <= p1;
      }

    // 2D cell
    else if (tdim == 2 && gdim == 2)
      {
        /*v1:
        auto v0 = v(0);
        auto v1 = v(1);
        auto v2 = v(2);
        */
        auto v0 = vertices[0];
        auto v1 = vertices[1];
        auto v2 = vertices[2];

        // mesh.coordinates() returns a 1D array of doubles: x0, y0, x1, y1,... in
        // this case.  So the value of x for vertex n is at offset n*2. Set pointers
        // to each vertex point:
        double* p0 = &mesh.coordinates()[v0*gdim];
        double* p1 = &mesh.coordinates()[v1*gdim];
        double* p2 = &mesh.coordinates()[v2*gdim];

        // The rest is code copied from CollisionPredicates.cpp
  
        const double ref = dnt::_orient2d(p0, p1, p2);

        if (ref > 0.0)
          {
            return (dnt::_orient2d(p1, p2, point) >= 0.0 and
                    dnt::_orient2d(p2, p0, point) >= 0.0 and
                    dnt::_orient2d(p0, p1, point) >= 0.0);
          }
        else if (ref < 0.0)
          {
            return (dnt::_orient2d(p1, p2, point) <= 0.0 and
                    dnt::_orient2d(p2, p0, point) <= 0.0 and
                    dnt::_orient2d(p0, p1, point) <= 0.0);
          }
        else
          {
            return ((dnt::_orient2d(p0, p1, point) == 0.0 and
                     dnt::_collides_segment_point_1d(p0[0], p1[0], point[0]) and
                     dnt::_collides_segment_point_1d(p0[1], p1[1], point[1])) or
                    (dnt::_orient2d(p1, p2, point) == 0.0 and
                     dnt::_collides_segment_point_1d(p1[0], p2[0], point[0]) and
                     dnt::_collides_segment_point_1d(p1[1], p2[1], point[1])) or
                    (dnt::_orient2d(p2, p0, point) == 0.0 and
                     dnt::_collides_segment_point_1d(p2[0], p0[0], point[0]) and
                     dnt::_collides_segment_point_1d(p2[1], p0[1], point[1])));
          }
      }

    // 3D cell
    else if (tdim == 3 && gdim == 3)
      {
        auto v0 = vertices[0];
        auto v1 = vertices[1];
        auto v2 = vertices[2];
        auto v3 = vertices[3];

        //  std::cout << "v0: " << v0 << " v1: " << v1 << " v2: " << v2 << " v3: " << v3 << std::endl;
        //  std::cout << "x: " << x << " y: " << y << " z: " << z << std::endl;
  
        // mesh.coordinates() returns a 1D array of doubles: x0, y0, z0, x1, y1, z1,... in this
        // case.  So the value of x for vertex n is at n*2.
        double* p0 = &mesh.coordinates()[v0*gdim];
        double* p1 = &mesh.coordinates()[v1*gdim];
        double* p2 = &mesh.coordinates()[v2*gdim];
        double* p3 = &mesh.coordinates()[v3*gdim];

        // std::cout << "p0: " << p0[0] << ", " << p0[1] << ", " << p0[2] << "\n"
        //           << "p1: " << p1[0] << ", " << p1[1] << ", " << p1[2] << "\n"
        //           << "p2: " << p2[0] << ", " << p2[1] << ", " << p2[2] << "\n"
        //           << "p3: " << p3[0] << ", " << p3[1] << ", " << p3[2]
        //           << std::endl;

        //  std::cout << "Point coords p: " << point[0] << " py: " << point[1] <<" pz: " << point[2] << std::endl;

        const double ref = dnt::_orient3d(p0, p1, p2, p3);

        if (ref > 0.0)
          {
            return (dnt::_orient3d(p0, p1, p2, point) >= 0.0 and
                    dnt::_orient3d(p0, p3, p1, point) >= 0.0 and
                    dnt::_orient3d(p0, p2, p3, point) >= 0.0 and
                    dnt::_orient3d(p1, p3, p2, point) >= 0.0);
          }
        else if (ref < 0.0)
          {
            return (dnt::_orient3d(p0, p1, p2, point) <= 0.0 and
                    dnt::_orient3d(p0, p3, p1, point) <= 0.0 and
                    dnt::_orient3d(p0, p2, p3, point) <= 0.0 and
                    dnt::_orient3d(p1, p3, p2, point) <= 0.0);
          }
        else
          {
            dolfin::dolfin_error("cell_contains_point (3D)",
                                 "compute tetrahedron point collision",
                                 "Not implemented for degenerate tetrahedron");
          }
      }

    return false;
  }    
  

  // We moved these to dolfin_cpp.cpp

  /*
    void dolfin_DnT_pstruct_instances()
    {
    dolfin::Mesh mesh;
    dolfin::Mesh& meshRef = mesh;
    py::array_t<int> vertices;
    Pstruct<Ptype::cartesian_x> ps1D;
    Pstruct<Ptype::cartesian_xy> ps2D;
    Pstruct<Ptype::cartesian_xyz> ps3D;

    // c++filt said that the following symbols were not in the .so file.  These
    // statements force the compiler to make these specialized versions of the
    // templated function in dolfin.o.

    // NEEDED?  
    cell_contains_point(meshRef, vertices, ps1D);
    cell_contains_point(meshRef, vertices, ps2D);
    cell_contains_point(meshRef, vertices, ps3D);
    }
  */

  //! Find the cell facet crossed in traveling along a displacement vector dr from position r0

  /*!
    Overload on the number of particle coordinates: 1D, 2D, 3D 
  */

  /// Versions using DnT_struct

  // 1D version using a DnT_pstruct
  //bool cell_contains_point(dolfin::Mesh& mesh,
  //                         py::array_t<int> vertices,
  //                         Pstruct<Ptype::cartesian_x> ps1D)

} // namespace dnt
