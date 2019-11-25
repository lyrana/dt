/*! \file dolfin_functions.cpp

  \brief This file has the source code for a C++ implementation of dolfin-based functions.

  The following C++ functions are defined:
    divide_by_cell_volumes()
    interpolate_field_to_points()
    is_inside_vertices()
    find_facet()
    XXparticle_stays_in_cell()XX

  \namespace dnt
  \sa dolfin_functions.h, dolfin_functions_cpp.cpp

*/

#include "MeshEntityArrays.h"
#include "dolfin_functions.h"
#include "predicates.h"

// uncomment to disable assert()
// #define NDEBUG
#include <cassert>

namespace py = pybind11;

/*
void init_ex2(py::module &m) {
m.def("sub", [](int a, int b) { return a - b; });
}
*/

// Static scratch space local to this file.
std::vector<double> vertex_coords(12);
double vecToFacet[3];

namespace dnt
{
  // BEGINDEF: bool is_inside_vertices(dolfin::Mesh& mesh, const unsigned int* vertices, double* point)
  //! Test whether a point lies within the cell defined by the vertices.
  /*!

  */
  bool is_inside_vertices(dolfin::Mesh& mesh,
                           const unsigned int* vertices,
                           double* point)
  {
    // The algorithm used depends on the mesh dimension
    const std::size_t tDim = mesh.topology().dim();
    const std::size_t gDim = mesh.geometry().dim();

    //v1:    auto v = vertices.unchecked<1>(); // <1> is the number of array indices. See pybind11 docs.

    // 1D mesh
    if (tDim == 1 && gDim == 1)
      {
        //v1:        auto v0 = v(0);
        //v1:        auto v1 = v(1);
        auto v0 = vertices[0];
        auto v1 = vertices[1];
        double p0 = mesh.coordinates()[v0];
        double p1 = mesh.coordinates()[v1];

        //        std::cout << "dolfin_functions.cpp::is_inside_vertices: point[0] = " << point[0] << " vertex p0 = " << p0 << ", p1 = " << p1 << std::endl;

        if (p0 > p1)
          std::swap(p0, p1);
        return p0 <= point[0] and point[0] <= p1;
      }

    // 2D mesh
    else if (tDim == 2 && gDim == 2)
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
        double* p0 = &mesh.coordinates()[v0*gDim];
        double* p1 = &mesh.coordinates()[v1*gDim];
        double* p2 = &mesh.coordinates()[v2*gDim];

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

    // 3D mesh
    else if (tDim == 3 && gDim == 3)
      {
        auto v0 = vertices[0];
        auto v1 = vertices[1];
        auto v2 = vertices[2];
        auto v3 = vertices[3];

        //  std::cout << "v0: " << v0 << " v1: " << v1 << " v2: " << v2 << " v3: " << v3 << std::endl;
        //  std::cout << "x: " << x << " y: " << y << " z: " << z << std::endl;
  
        // mesh.coordinates() returns a 1D array of doubles: x0, y0, z0, x1, y1, z1,... in this
        // case.  So the value of x for vertex n is at n*2.
        double* p0 = &mesh.coordinates()[v0*gDim];
        double* p1 = &mesh.coordinates()[v1*gDim];
        double* p2 = &mesh.coordinates()[v2*gDim];
        double* p3 = &mesh.coordinates()[v3*gDim];

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
            dolfin::dolfin_error("is_inside_vertices (3D)",
                                 "compute tetrahedron point collision",
                                 "Not implemented for degenerate tetrahedron");
          }
      }

    return false;
  } // ENDDEF: bool is_inside_vertices(dolfin::Mesh& mesh, const unsigned int* vertices, double* point)

  
  // BEGINDEF: py::tuple find_facet(py::object mesh_M, double* x0, double* dx, size_t cell_index)
  //! Find the cell facet crossed in traveling along a displacement vector dr from position r0.
  /*!
    The facet crossed is the one with the smallest value of the
    fraction (normal distance to facet plane)/(total distance
    traveled normal to facet plane)

    1. Compute the distance a particle moves in the directions
       of a facet (Lf).

    2. If it moved closer to it, compute the distance to the
       facet (Df).

    3. If Df < Lf, the particle crossed the plane that the
       facet is in.

    4. If the particle crosses more than one facet plane, the
       facet crossed is the one with the smallest ratio Df/Lf,
       since it reached that facet first.

       \param x0: initial position, as a tuple of coordinates (x),
                  (x,y), or (x, y, z)

       \param dx: move vector, as a tuple of coordinates (x),
                  (x,y), or (x, y, z)
           
       \param cell_index: the (global?) index of the cell that the
                          particle started in.

       \returns: 3-tuple(index of the facet crossed or None,
                         fraction of the path that's in this cell,
                         the unit normal to the facet crossed)
  */
  py::tuple find_facet(py::object mesh_M, double* x0, double* dx, size_t cell_index)
  {
    // Get attributes needed from Mesh_C arg
    auto mesh = mesh_M.attr("mesh").cast<dolfin::Mesh>();
    /*
    auto cellFacetNormalsDictCap = mesh_M.attr("cell_facet_normals_dict").cast<py::capsule>();
    std::cout << "find_facet(): The name of the capsule is " << cellFacetNormalsDictCap.name() << std::endl;        
    std::map<int, double*>* cellFacetNormalsDict = cellFacetNormalsDictCap;
    */
    // The facet-normals are in mesh_M.attr("mea_object"). The array
    // containing them is templated on the number of facets in the cell, so they are
    // accessed below, depending on the cell type. The vectors are laid out as n0,
    // n1, ... to the number of facets, where the n's are the normals. The normals
    // are 3-vectors, regardless of the geometric dimension.
    size_t cellFNVdim = 3; // Number of doubles per facet-normal vector (FNV).
    
    auto NO_FACET = mesh_M.attr("NO_FACET").cast<int>(); // NO_FACET is a static class constant
                                                         // of Mesh_C defined in Dolfin_Module.py
    const std::size_t gDim = mesh.geometry().dim(); // The geometric dimension of the mesh

    int facet(NO_FACET);
    double dxFraction(1.0);

    double dxDotdx(0.0);
    if (gDim == 3)
      {
        dxDotdx = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
      }
    else if (gDim == 2)
      {
        dxDotdx = dx[0]*dx[0] + dx[1]*dx[1];
      }
    else if (gDim == 1)
      {
        dxDotdx = dx[0]*dx[0];
      }
  
    if (dxDotdx == 0.0)
      {
        return py::make_tuple(facet, dxFraction, nullptr);
      }

    // cell = self.cell_dict[cell_index];
    dolfin::Cell cell(mesh, cell_index);
    
    // The coordinates of all the vertices in this cell.

    //Cell.h:339: void get_vertex_coordinates(std::vector<double>& coordinates) const
      // for (std::size_t i = 0; i < num_vertices; i++)
      //   for (std::size_t j = 0; j < gDim; j++)
      //     coordinates[i*gDim + j] = _mesh->geometry().x(vertices[i])[j];
    // The coordinates are laid out as (x0,y0,z0), (x1,y1,z1), ... for the vertices
    cell.get_vertex_coordinates(vertex_coords);
    
    if (gDim == 3)
      {
        // 3D mesh: There are 4 facets indexed 0,1,2,3.

        // The facet-normals are laid out as n0, n1, n2, n3, where
        // the n's are the normals. The normals are 3-vectors.
        // The type here is std::array<double, N_CELL_FACETS*3>. The 3 is for 3-vector.
        auto meshEntityArrays = mesh_M.attr("mea_object").cast<MeshEntityArrays<4> *>();
        auto facetNormalVectors = meshEntityArrays->get_cell_facet_normals(cell_index);

        // Test if the plane of facet 0 is crossed:
        auto n0 = &facetNormalVectors[0];
        // n0DotDx = np_m.dot(facetNormalVectors[0], dx)
        auto n0DotDx = n0[0]*dx[0] + n0[1]*dx[1] + n0[2]*dx[2];
        
        if (n0DotDx > 0)
          {
            // The vector from the starting point to a vertex in the facet plane. facet 0 is opposite vertex 0, so use vertex 1 here:
            // vecToFacet = np_m.subtract(vertex_coords[1], x0)
            auto vertex1 = &vertex_coords[1*gDim]; // Coordinates for vertex 1
            vecToFacet[0] = vertex1[0] - x0[0];
            vecToFacet[1] = vertex1[1] - x0[1];
            vecToFacet[2] = vertex1[2] - x0[2];
              
            // The normal distance to the facet plane
            // distanceToFacet = np_m.dot(facetNormalVectors[0], vecToFacet)
            auto distanceToFacet = n0[0]*vecToFacet[0] + n0[1]*vecToFacet[1] + n0[2]*vecToFacet[2];
//                print "f0 find_facet(): vecToFacet=", vecToFacet, "distanceToFacet=", distanceToFacet
            if (distanceToFacet < 0.0) // Assume this is due to round-off error and flip the sign
              {
                std::cout << "dolfin_functions.cpp::find_facet 0: !!! Bad value for distanceToFacet: " << distanceToFacet << ". Assuming it's a tiny number and flipping the sign to continue!!!" << std::endl;
              distanceToFacet = -distanceToFacet;
              }
            if (distanceToFacet < dxFraction*n0DotDx)
              {
                facet = 0; // The plane of facet 0 was crossed
                dxFraction = distanceToFacet/n0DotDx;
              }
          }

        // Next, test if the plane of facet 1 is crossed:
        auto n1 = &facetNormalVectors[1*cellFNVdim];
        auto n1DotDx = n1[0]*dx[0] + n1[1]*dx[1] + n1[2]*dx[2];
        if (n1DotDx > 0)
          {
            auto vertex2 = &vertex_coords[2*gDim]; // Coordinates for vertex 2
            vecToFacet[0] = vertex2[0] - x0[0];
            vecToFacet[1] = vertex2[1] - x0[1];
            vecToFacet[2] = vertex2[2] - x0[2];
              
            // The normal distance to the facet plane
            // distanceToFacet = np_m.dot(facetNormalVectors[0], vecToFacet)
            auto distanceToFacet = n1[0]*vecToFacet[0] + n1[1]*vecToFacet[1] + n1[2]*vecToFacet[2];
//                print "f1 find_facet(): vecToFacet=", vecToFacet, "distanceToFacet=", distanceToFacet
            if (distanceToFacet < 0.0) // Assume this is due to round-off error and flip the sign
              {
                std::cout << "dolfin_functions.cpp::find_facet 1: !!! Bad value for distanceToFacet: " << distanceToFacet << ". Assuming it's a tiny number and flipping the sign to continue!!!" << std::endl;
              distanceToFacet = -distanceToFacet;
              }
            if (distanceToFacet < dxFraction*n1DotDx)
              {
                facet = 1; // The plane of facet 1 was crossed
                dxFraction = distanceToFacet/n1DotDx;
              }
          }

        // Next, test if the plane of facet 2 is crossed:
        auto n2 = &facetNormalVectors[2*cellFNVdim];
        auto n2DotDx = n2[0]*dx[0] + n2[1]*dx[1] + n2[2]*dx[2];
        if (n2DotDx > 0)
          {
            auto vertex3 = &vertex_coords[3*gDim]; // Coordinates for vertex 3
            vecToFacet[0] = vertex3[0] - x0[0];
            vecToFacet[1] = vertex3[1] - x0[1];
            vecToFacet[2] = vertex3[2] - x0[2];
              
            // The normal distance to the facet plane
            // distanceToFacet = np_m.dot(facetNormalVectors[0], vecToFacet)
            auto distanceToFacet = n2[0]*vecToFacet[0] + n2[1]*vecToFacet[1] + n2[2]*vecToFacet[2];
//                print "f2 find_facet(): vecToFacet=", vecToFacet, "distanceToFacet=", distanceToFacet
            if (distanceToFacet < 0.0) // Assume this is due to round-off error and flip the sign
              {
                std::cout << "dolfin_functions.cpp::find_facet 2: !!! Bad value for distanceToFacet: " << distanceToFacet << ". Assuming it's a tiny number and flipping the sign to continue!!!" << std::endl;
              distanceToFacet = -distanceToFacet;
              }
            if (distanceToFacet < dxFraction*n2DotDx)
              {
                facet = 2; // The plane of facet 2 was crossed
                dxFraction = distanceToFacet/n2DotDx;
              }
          }
        
        // Next, test if the plane of facet 3 is crossed:
        auto n3 = &facetNormalVectors[3*cellFNVdim];
        auto n3DotDx = n3[0]*dx[0] + n3[1]*dx[1] + n3[2]*dx[2];
        if (n3DotDx > 0)
          {
            auto vertex0 = &vertex_coords[0]; // Coordinates for vertex 0 (opposite facet 3)
            vecToFacet[0] = vertex0[0] - x0[0];
            vecToFacet[1] = vertex0[1] - x0[1];
            vecToFacet[2] = vertex0[2] - x0[2];
              
            // The normal distance to the facet plane
            // distanceToFacet = np_m.dot(facetNormalVectors[0], vecToFacet)
            auto distanceToFacet = n3[0]*vecToFacet[0] + n3[1]*vecToFacet[1] + n3[2]*vecToFacet[2];
//                print "f3 find_facet(): vecToFacet=", vecToFacet, "distanceToFacet=", distanceToFacet
            if (distanceToFacet < 0.0) // Assume this is due to round-off error and flip the sign
              {
                std::cout << "dolfin_functions.cpp::find_facet 3: !!! Bad value for distanceToFacet: " << distanceToFacet << ". Assuming it's a tiny number and flipping the sign to continue!!!" << std::endl;
              distanceToFacet = -distanceToFacet;
              }
            if (distanceToFacet < dxFraction*n3DotDx)
              {
                facet = 3; // The plane of facet 3 was crossed
                dxFraction = distanceToFacet/n3DotDx;
              }
          }
        // Check range of dxFraction:
        if (dxFraction > 1.0 || dxFraction < 0.0)
          {
            std::cout << "dolfin_functions.cpp::find_facet: !!! Bad value for dxFraction: " << dxFraction << std::endl;
            exit(EXIT_FAILURE);
          }
        
        return py::make_tuple(facet, dxFraction, facetNormalVectors[facet*cellFNVdim]);
                                     
      } // if (gDim == 3)
    else if (gDim == 2)
      {
        // 2D mesh: There are 3 facets indexed 0,1,2.

        // The facet-normals are laid out as n0, n1, n2, where
        // the n's are the normals. The normals are always 3-vectors.
        // The type here is std::array<double, N_CELL_FACETS*3>
        auto meshEntityArrays = mesh_M.attr("mea_object").cast<MeshEntityArrays<3> *>();
        auto facetNormalVectors = meshEntityArrays->get_cell_facet_normals(cell_index);
        
        // Test if the plane of facet 0 is crossed:
        auto n0 = &facetNormalVectors[0];
        auto n0DotDx = n0[0]*dx[0] + n0[1]*dx[1];
        
        if (n0DotDx > 0)
          {
            // The vector from the starting point to a vertex in the facet plane. facet 0 is opposite vertex 0, so use vertex 1 here:
            // vecToFacet = np_m.subtract(vertex_coords[1], x0)
            auto vertex1 = &vertex_coords[1*gDim]; // Coordinates for vertex 1
            vecToFacet[0] = vertex1[0] - x0[0];
            vecToFacet[1] = vertex1[1] - x0[1];
              
            // The normal distance to the facet plane
            // distanceToFacet = np_m.dot(facetNormalVectors[0], vecToFacet)
            auto distanceToFacet = n0[0]*vecToFacet[0] + n0[1]*vecToFacet[1];
//                print "f0 find_facet(): vecToFacet=", vecToFacet, "distanceToFacet=", distanceToFacet
            if (distanceToFacet < 0.0) // Assume this is due to round-off error and flip the sign
              {
                std::cout << "dolfin_functions.cpp::find_facet 0: !!! Bad value for distanceToFacet: " << distanceToFacet << ". Assuming it's a tiny number and flipping the sign to continue!!!" << std::endl;
              distanceToFacet = -distanceToFacet;
              }
            if (distanceToFacet < dxFraction*n0DotDx)
              {
                facet = 0; // The plane of facet 0 was crossed
                dxFraction = distanceToFacet/n0DotDx;
              }
          }

        // Next, test if the plane of facet 1 is crossed:
        auto n1 = &facetNormalVectors[1*cellFNVdim];
        auto n1DotDx = n1[0]*dx[0] + n1[1]*dx[1];
        if (n1DotDx > 0)
          {
            auto vertex2 = &vertex_coords[2*gDim]; // Coordinates for vertex 2
            vecToFacet[0] = vertex2[0] - x0[0];
            vecToFacet[1] = vertex2[1] - x0[1];
              
            // The normal distance to the facet plane
            // distanceToFacet = np_m.dot(facetNormalVectors[0], vecToFacet)
            auto distanceToFacet = n1[0]*vecToFacet[0] + n1[1]*vecToFacet[1];
//                print "f1 find_facet(): vecToFacet=", vecToFacet, "distanceToFacet=", distanceToFacet
            if (distanceToFacet < 0.0) // Assume this is due to round-off error and flip the sign
              {
                std::cout << "dolfin_functions.cpp::find_facet 1: !!! Bad value for distanceToFacet: " << distanceToFacet << ". Assuming it's a tiny number and flipping the sign to continue!!!" << std::endl;
              distanceToFacet = -distanceToFacet;
              }
            if (distanceToFacet < dxFraction*n1DotDx)
              {
                facet = 1; // The plane of facet 1 was crossed
                dxFraction = distanceToFacet/n1DotDx;
              }
          }

        // Next, test if the plane of facet 2 is crossed:
        auto n2 = &facetNormalVectors[2*cellFNVdim];
        auto n2DotDx = n2[0]*dx[0] + n2[1]*dx[1];
        if (n2DotDx > 0)
          {
            auto vertex0 = &vertex_coords[0]; // Coordinates for vertex 0 (opposite facet 2)
            vecToFacet[0] = vertex0[0] - x0[0];
            vecToFacet[1] = vertex0[1] - x0[1];
              
            // The normal distance to the facet plane
            // distanceToFacet = np_m.dot(facetNormalVectors[0], vecToFacet)
            auto distanceToFacet = n2[0]*vecToFacet[0] + n2[1]*vecToFacet[1];
//                print "f2 find_facet(): vecToFacet=", vecToFacet, "distanceToFacet=", distanceToFacet
            if (distanceToFacet < 0.0) // Assume this is due to round-off error and flip the sign
              {
                std::cout << "dolfin_functions.cpp::find_facet 2: !!! Bad value for distanceToFacet: " << distanceToFacet << ". Assuming it's a tiny number and flipping the sign to continue!!!" << std::endl;
              distanceToFacet = -distanceToFacet;
              }
            if (distanceToFacet < dxFraction*n2DotDx)
              {
                facet = 2; // The plane of facet 2 was crossed
                dxFraction = distanceToFacet/n2DotDx;
              }
          }
        
        // Check range of dxFraction:
        if (dxFraction > 1.0 || dxFraction < 0.0)
          {
            std::cout << "dolfin_functions.cpp::find_facet: !!! Bad value for dxFraction: " << dxFraction << std::endl;
            exit(EXIT_FAILURE);
          }
        
        return py::make_tuple(facet, dxFraction, facetNormalVectors[facet*cellFNVdim]);
                                     
      } // else if (gDim == 2)
    else if (gDim == 1)
      {
        // 1D mesh: There are 2 facets, indexed 0,1.
        // In 1D, facet 0 and node 0 have the same index, 0.

        // The facet-normals are laid out as n0, n1 where
        // the n's are the normals. The normals are always 3-vectors.
        // The type here is std::array<double, N_CELL_FACETS*3>
        auto meshEntityArrays = mesh_M.attr("mea_object").cast<MeshEntityArrays<2> *>();
        auto facetNormalVectors = meshEntityArrays->get_cell_facet_normals(cell_index);

        // Test if the plane of facet 0 is crossed:
        auto n0 = &facetNormalVectors[0];
        auto n0DotDx = n0[0]*dx[0];

        //        std::cout << "find_facet: x0[0] " << x0[0] << " dx[0] " << dx[0] << " n0DotDx " << n0DotDx << std::endl;
        
        if (n0DotDx > 0)
          {
            // The vector from the starting point to a vertex in the facet plane.
            auto vertex0 = &vertex_coords[0]; // Coordinates for vertex 0.
            vecToFacet[0] = vertex0[0] - x0[0];
               
            // The normal distance to the facet plane
            // distanceToFacet = np_m.dot(facetNormalVectors[0], vecToFacet)
            auto distanceToFacet = n0[0]*vecToFacet[0];
//                print "f0 find_facet(): vecToFacet=", vecToFacet, "distanceToFacet=", distanceToFacet
            if (distanceToFacet < 0.0) // Assume this is due to round-off error and flip the sign
              {
                std::cout << "dolfin_functions.cpp::find_facet 0: !!! Bad value for distanceToFacet: " << distanceToFacet << ". Assuming it's a tiny number and flipping the sign to continue!!!" << std::endl;
              distanceToFacet = -distanceToFacet;
              }
            if (distanceToFacet < dxFraction*n0DotDx)
              {
                facet = 0; // The plane of facet 0 was crossed
                dxFraction = distanceToFacet/n0DotDx;
                assert(dxFraction >= 0 && dxFraction <= 1.0);
                return py::make_tuple(facet, dxFraction, facetNormalVectors[facet*cellFNVdim]);
              }
          }

        // Next, test if the plane of facet 1 is crossed:
        auto n1 = &facetNormalVectors[1*cellFNVdim];
        auto n1DotDx = n1[0]*dx[0];
        if (n1DotDx > 0)
          {
            auto vertex1 = &vertex_coords[1*gDim]; // Coordinates for vertex 1
            vecToFacet[0] = vertex1[0] - x0[0];
              
            // The normal distance to the facet plane
            auto distanceToFacet = n1[0]*vecToFacet[0];
            //            std::cout << "f1 find_facet(): vecToFacet= " << vecToFacet[0] << " distanceToFacet= " << distanceToFacet << std::endl;
            if (distanceToFacet < 0.0) // Assume this is due to round-off error and flip the sign
              {
                std::cout << "dolfin_functions.cpp::find_facet 1: !!! Bad value for distanceToFacet: " << distanceToFacet << ". Assuming it's a tiny number and flipping the sign to continue!!!" << std::endl;
              distanceToFacet = -distanceToFacet;
              }
            if (distanceToFacet < dxFraction*n1DotDx)
              {
                facet = 1; // The plane of facet 1 was crossed
                dxFraction = distanceToFacet/n1DotDx;
                assert(dxFraction >= 0 && dxFraction <= 1.0);                
                return py::make_tuple(facet, dxFraction, facetNormalVectors[facet*cellFNVdim]);
              }
          }

        std::cout << "dolfin_functions.cpp::find_facet() for gDim = 1: Should never get here!!!" << std::endl;
        exit(EXIT_FAILURE);
        // return py::make_tuple(facet, dxFraction, facetNormalVectors[facet*cellFNVdim]);
                                     
      } // else if (gDim == 1)

    std::cout << "End of dolfin_functions.cpp::find_facet(). Should never get here!!!" << std::endl;
    exit(EXIT_FAILURE);
    
    return py::make_tuple(NO_FACET, 0.0, nullptr); // Shouldn't reach this. This stops a compiler warning.
    
  } // ENDDEF: py::tuple find_facet(py::object mesh_M, double* x0, double* dx, size_t cell_index)
  
  // We moved these to dolfin_functions_cpp.cpp

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
    // templated function in dolfin_functions.o.

    // NEEDED?  
    is_inside_vertices(meshRef, vertices, ps1D);
    is_inside_vertices(meshRef, vertices, ps2D);
    is_inside_vertices(meshRef, vertices, ps3D);
    }
  */

  //! Find the cell facet crossed in traveling along a displacement vector dr from position r0

  /*!
    Overload on the number of particle coordinates: 1D, 2D, 3D 
  */

  /// Versions using DnT_struct

  // 1D version using a DnT_pstruct
  //bool is_inside_vertices(dolfin::Mesh& mesh,
  //                         py::array_t<int> vertices,
  //                         Pstruct<Ptype::cartesian_x> ps1D)

} // namespace dnt
