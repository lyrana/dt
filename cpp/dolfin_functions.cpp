/*! \file dolfin_functions.cpp

  \brief This file has the source code for a C++ implementation of dolfin-based functions.

  The following C++ functions are defined:
    divide_by_cell_volumes()
    find_facet()
    interpolate_field_to_points()
      For py::array_t points arg
      For Pstruct* points arg
    is_inside_vertices()

  The Python bindings for these are in dolfin_functions_solib.cpp.

  The statements at the bottom of this file cause the compilation of the specialized
  versions needed in dolfin_functions_solib.cpp. Without these statements, those do not
  get compiled.

  \namespace dnt
  \sa dolfin_functions.h, MeshEntityArrays.h, predicates.h, dolfin_functions_solib.cpp

*/

#include "MeshEntityArrays.h"
#include "dolfin_functions.h"
#include "predicates.h"

// uncomment to disable assert()
// #define NDEBUG
#include <cassert>

namespace py = pybind11;

// Static scratch space local to this file.
std::vector<double> vertex_coords(12);
double vecToFacet[3];

namespace dnt
{
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

       \param returnStdArray: flag indicating whether the facet-normal is
                             returned as a Numpy array or as a pointer-to-data.

       \returns: 3-tuple(index of the facet crossed or None,
                         fraction of the path that's in this cell,
                         the unit normal to the facet crossed)
  */
  py::tuple find_facet(py::object mesh_M, double* x0, double* dx, size_t cell_index, bool returnStdArray)
  {

    //    std::cout << "Entering dolfin_functions.cpp::find_facet()." << std::endl;
    
    // Get attributes needed from Mesh_C arg
    auto mesh = mesh_M.attr("mesh").cast<dolfin::Mesh>();

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
    
    // Get the coordinates of all the vertices in this cell.
    // get_vertex_coordinates COPIES the vertex coordinates to the static array
    // vertex_coords[].

    //Cell.h:339: void get_vertex_coordinates(std::vector<double>& coordinates) const
      // for (std::size_t i = 0; i < num_vertices; i++)
      //   for (std::size_t j = 0; j < gDim; j++)
      //     coordinates[i*gDim + j] = _mesh->geometry().x(vertices[i])[j];
    // The coordinates are laid out as (x0,y0,z0), (x1,y1,z1), ... for the vertices
    cell.get_vertex_coordinates(vertex_coords);

    //    std::cout << "dolfin_functions.cpp::find_facet(). gDim " << gDim << std::endl;
    
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
        
        if (returnStdArray == true)
          {
            return py::make_tuple(facet, dxFraction, &facetNormalVectors[facet*cellFNVdim]);
          }
        else
          {
            // Create a Numpy array for the facet-normal and give it gDim values
            auto pyFacetNormalVector = py::array_t<double>(gDim, &facetNormalVectors[facet*cellFNVdim]);
            return py::make_tuple(facet, dxFraction, pyFacetNormalVector);
          }
                                     
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
        if (returnStdArray == true)
          {
            return py::make_tuple(facet, dxFraction, &facetNormalVectors[facet*cellFNVdim]);
          }
        else
          {
            // Create a Numpy array for the facet-normal and give it gDim values
            auto pyFacetNormalVector = py::array_t<double>(gDim, &facetNormalVectors[facet*cellFNVdim]);
            return py::make_tuple(facet, dxFraction, pyFacetNormalVector);
          }
                                     
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

        //        std::cout << "find_facet (gDim = 1): x0[0] " << x0[0] << " dx[0] " << dx[0] << " n0DotDx " << n0DotDx << std::endl;
        
        if (n0DotDx > 0)
          {
            // The vector from the starting point to a vertex in the facet plane.
            auto vertex0 = &vertex_coords[0]; // Coordinates for vertex 0.
            vecToFacet[0] = vertex0[0] - x0[0];
               
            // The normal distance to the facet plane
            // distanceToFacet = np_m.dot(facetNormalVectors[0], vecToFacet)
            auto distanceToFacet = n0[0]*vecToFacet[0];
            //            std::cout << "f0 find_facet(): vecToFacet= " << vecToFacet[0] <<  " distanceToFacet= " << distanceToFacet << std::endl;
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
                //                std::cout << "dolfin_functions.cpp::find_facet0. facetNormalVectors[facet*cellFNVdim] is " << facetNormalVectors[facet*cellFNVdim] << std::endl;                

                if (returnStdArray == true)
                  {
                    return py::make_tuple(facet, dxFraction, &facetNormalVectors[facet*cellFNVdim]);
                  }
                else
                  {
                    // Create a Numpy array for the facet-normal and give it gDim values                    
                    auto pyFacetNormalVector = py::array_t<double>(gDim, &facetNormalVectors[facet*cellFNVdim]);
                    return py::make_tuple(facet, dxFraction, pyFacetNormalVector);
                  }
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
                if (returnStdArray == true)
                  {
                    return py::make_tuple(facet, dxFraction, &facetNormalVectors[facet*cellFNVdim]);
                  }
                else
                  {
                    // Create a Numpy array for the facet-normal and give it gDim values
                    auto pyFacetNormalVector = py::array_t<double>(gDim, &facetNormalVectors[facet*cellFNVdim]);
                    return py::make_tuple(facet, dxFraction, pyFacetNormalVector);
                  }
              }
          }

        std::cout << "dolfin_functions.cpp::find_facet() for gDim = 1: Should never get here!!!" << std::endl;
        exit(EXIT_FAILURE);
                                     
      } // else if (gDim == 1)

    std::cout << "End of dolfin_functions.cpp::find_facet(). Should never get here!!!" << std::endl;
    exit(EXIT_FAILURE);
    
    return py::make_tuple(NO_FACET, 0.0, nullptr); // Shouldn't reach this. This statement just stops a compiler warning.
    
  } // ENDDEF: py::tuple find_facet(py::object mesh_M, double* x0, double* dx, size_t cell_index, bool returnStdArray)


  //! (Version for py::array_t points) Evaluate a dolfin Function at a set of points. The point data are in a structured Numpy array.
  /*!

    This is the C++ version of the Python function in Dolfin_Module.py:
    interpolate_field_to_points(self, points, field_at_points).

    Note that, because E is constant-over-element, no linear interpolation is needed.

    \param field is a reference to a Dolfin Function that has the values of the
           field on a finite-element mesh.

    \param points is a structured Numpy array containing the point data.

    \param npoints is the number of points in "points" to interpolate to.

    \param field_at_points is a structured Numpy array to hold the values computed at the
           points. This is laid out in C-order as [number_of_points, number_of_field_components].

     The "points" Numpy array is accessed using a "unchecked proxy" to bypass
     array bounds-checking.

  */
  //  template <Ptype PT, typename Ftype>
  template <Ptype PT>
  void interpolate_field_to_points(dolfin::Function* field,
                                   py::array_t<Pstruct<PT>, 0> points,
                                   py::ssize_t npoints,
                                   py::array_t<double> field_at_points)
  {
    //std::cout << "Hello from interpolate_field_to_points()@A1"  << std::endl;
    
    // Use the buffer protocol to access the data in points:
    /*
    const auto points_info = points.request(); // request() returns metadata about the Python array (ptr, ndim, size, shape)
    const auto pArray = static_cast<Pstruct<PT>*>(points_info.ptr); // Pointer to a struct in points
    */
    // Use an unchecked proxy instead:
    auto pointsProxy = points.unchecked();
    
    std::shared_ptr< const dolfin::FunctionSpace > functionSpace = field->function_space();

    // The value of dofs_per_cell is the number of E-components times the order-of-the-element (1 for constant E-in-element)
    // std::size_t dofs_per_cell = functionSpace->element()->space_dimension();
    // std::cout << "interpolate_field_to_points(): dofs_per_cell = " << dofs_per_cell << std::endl;
    std::shared_ptr< const dolfin::GenericDofMap > dofMap = functionSpace->dofmap();
  
    auto fieldAtPointsProxy = field_at_points.mutable_unchecked<2>(); // fieldAtPointsProxy is a 2D array of doubles.
    // std::cout << " shape(0)= " << fieldAtPointsProxy.shape(0) << " shape(1)= " << fieldAtPointsProxy.shape(1) << std::endl;

    // Pull out the components of the (constant) field in this cell
    // auto nComps = fieldAtPointsProxy.shape(1);
    for (auto ip = 0; ip < npoints; ip++) {
      // std::cout << "particle " << ip << " x= " << pointsProxy[ip].x_ << std::endl;
      
      auto cellIndex = pointsProxy(ip).cell_index_;
      
      // std::cout << "particle " << ip << " cellIndex= " << cellIndex << std::endl;
      // In the case of constant-E-in-cell, there's just one set of E values in the cell, so we
      // use the name 'fieldValue' instead of 'fieldValues':
      auto fieldValue = dofMap->cell_dofs(cellIndex);
      // std::cout << "particle " << ip << " nComps= " << fieldValue.size() << std::endl;
      field->vector()->get_local(&fieldAtPointsProxy(ip, 0), fieldValue.size(), fieldValue.data());
      
      // Copy the field vector values to the field struct
      // double_to_fstruct<Ftype>(values, field_at_points[ip]);
    }

  }
  //ENDDEF: void interpolate_field_to_points(dolfin::Function& field,...

    //! (Version for Pstruct* points) Evaluate a dolfin Function at a set of points. The point data are in a structured Numpy array.
  /*!

    This is the C++ version of the Python function in Dolfin_Module.py:
    interpolate_field_to_points(self, points, field_at_points).

    Note that, because E is constant-over-element, no linear interpolation is needed.

    \param field is a reference to a Dolfin Function that has the values of the
           field on a finite-element mesh.

    \param points is a standard C array of Pstructs containing the point data.

    \param npoints is the number of points in "points" to interpolate to.

    \param field_at_points is a structured Numpy array to hold the values computed at the
           points. This is laid out in C-order as [number_of_points, number_of_field_components].

  */
  //  template <Ptype PT, typename Ftype>
  template <Ptype PT>
  void interpolate_field_to_points(dolfin::Function* field,
                                   Pstruct<PT>* points,
                                   py::ssize_t npoints,
                                   py::array_t<double> field_at_points)
  {
    // std::cout << "Hello from interpolate_field_to_points()@B1"  << std::endl;

    std::shared_ptr< const dolfin::FunctionSpace > functionSpace = field->function_space();

    // The value of dofs_per_cell is the number of E-components times the order-of-the-element (1 for constant E-in-element)
    // std::size_t dofs_per_cell = functionSpace->element()->space_dimension();
    // std::cout << "interpolate_field_to_points(): dofs_per_cell = " << dofs_per_cell << std::endl;
    std::shared_ptr< const dolfin::GenericDofMap > dofMap = functionSpace->dofmap();
  
    auto fieldAtPointsProxy = field_at_points.mutable_unchecked<2>(); // fieldAtPointsProxy is a 2D array of doubles.
    // std::cout << " shape(0)= " << fieldAtPointsProxy.shape(0) << " shape(1)= " << fieldAtPointsProxy.shape(1) << std::endl;

    // Pull out the components of the (constant) field in this cell
    // auto nComps = fieldAtPointsProxy.shape(1);
    for (auto ip = 0; ip < npoints; ip++) {
      // std::cout << "particle " << ip << " x= " << pointsProxy[ip].x_ << std::endl;
      
      auto cellIndex = points[ip].cell_index_;      
      
      // std::cout << "particle " << ip << " cellIndex= " << cellIndex << std::endl;
      // In the case of constant-E-in-cell, there's just one set of E values in the cell, so we
      // use the name 'fieldValue' instead of 'fieldValues':
      auto fieldValue = dofMap->cell_dofs(cellIndex);
      // std::cout << "particle " << ip << " nComps= " << fieldValue.size() << std::endl;
      field->vector()->get_local(&fieldAtPointsProxy(ip, 0), fieldValue.size(), fieldValue.data());
      
      // Copy the field vector values to the field struct
      // double_to_fstruct<Ftype>(values, field_at_points[ip]);
    }

  }
  //ENDDEF: void interpolate_field_to_points(dolfin::Function& field,...


  //! Test whether a point lies within the simplex (cell) defined by the vertices.
  /*!

  */
  bool is_inside_vertices(dolfin::Mesh& mesh,
                          const unsigned int* vertices,
                          double* point)
  {
    // The algorithm used depends on the mesh dimension
    const std::size_t tDim = mesh.topology().dim();
    const std::size_t gDim = mesh.geometry().dim();

    // 1D mesh
    if (tDim == 1 && gDim == 1)
      {
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

  
} // namespace dnt

// template void dnt::interpolate_field_to_points<>(dolfin::Function& field,
//                                                  py::array_t<Pstruct<Ptype::cartesian_xyz>, 0> points,
//                                                  py::ssize_t npoints,
//                                                  py::array_t<double> field_at_points);


// From https://en.cppreference.com/w/cpp/language/function_template:
// A function template by itself is not a type, or a function, or any other entity. No
// code is generated from a source file that contains only template definitions. In order
// for any code to appear, a template must be instantiated: the template arguments must be
// determined so that the compiler can generate an actual function (or class, from a class
// template).


// Since dolfin_functions_solib.cpp includes dolfin_function.h, but not dolfin_functions_cpp,
// we to need the following statements instantiating the function
// interpolate_field_to_points() with a particular type, to make the compiler compile the
// function.

// Use this one if given the points in a Numpy array
template void dnt::interpolate_field_to_points<>(dolfin::Function*,
                                                 py::array_t<Pstruct<Ptype::cartesian_xyz>, 0>,
                                                 py::ssize_t,
                                                 py::array_t<double>);
template void dnt::interpolate_field_to_points<>(dolfin::Function*,
                                                 py::array_t<Pstruct<Ptype::cartesian_xy>, 0>,
                                                 py::ssize_t,
                                                 py::array_t<double>);

// Use this one if given a pointer to the particle-struct array
template void dnt::interpolate_field_to_points<>(dolfin::Function*,
                                                 Pstruct<Ptype::cartesian_xyz>* points,
                                                 py::ssize_t,
                                                 py::array_t<double>);
template void dnt::interpolate_field_to_points<>(dolfin::Function*,
                                                 Pstruct<Ptype::cartesian_xy>* points,
                                                 py::ssize_t,
                                                 py::array_t<double>);
