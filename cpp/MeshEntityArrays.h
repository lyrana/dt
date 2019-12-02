/*! \file MeshEntityArrays.h

  \brief This is the header file for MeshEntityArrays.cpp

  \namespace dnt
  \sa MeshEntityArrays.cpp

*/

#ifndef MESHENTITYARRAYS_H
#define MESHENTITYARRAYS_H

#include <dolfin/common/Array.h>
#include <dolfin/common/Hierarchical.h>
#include <dolfin/function/assign.h>
#include <dolfin/function/Constant.h>
#include <dolfin/function/Expression.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionAssigner.h>
#include <dolfin/function/FunctionAXPY.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/MultiMeshFunction.h>
#include <dolfin/function/MultiMeshFunctionSpace.h>
#include <dolfin/function/MultiMeshSubSpace.h>
#include <dolfin/function/LagrangeInterpolator.h>
#include <dolfin/function/SpecialFunctions.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/fem/MultiMeshDofMap.h>
#include <dolfin/geometry/Point.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Facet.h>
//#include <pybind11/pybind11.h>

#include <Eigen/Dense>

/*
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
// #include <pybind11/eigen.h>
#include <pybind11/operators.h>
//#include <pybind11/cast.h>

namespace py = pybind11;
*/


// Put these in namespace dnt:

namespace dnt
{
  
/*! \class MeshEntityArrays
    \brief The MeshEntityArrays class provides C++ versions of mesh-entity dictionaries for moving particles.

    MeshEntityArrays is templated on the number of facets that a cell has. This is so
    that an std::array of appropriate size can be used to hold the facet-normal
    vectors.

    This class implements these mesh-entity dictionaries:

      cell_neighbors_array: This is a 2D boost:multi_array of ints.
      cell_neighbors_array[i][j] is the index of the cell attached to facet j
      of cell i. The int type is used because an index of -1 (NO_CELL) is used to
      indicate a boundary.

      cell_facet_normals_array: This is an std::vector containing the facet-normal
      vectors for each cell in the mesh, in an std::array<double, N_CELL_FACETS*3>.
      cell_facet_normals_array[i] contains the array [n0, n1, ...], where n0, n1,
      ... are the unit normal 3-vectors on each facet of cell i.
      The normal vector always has 3 components, regardless of the dimension of the mesh.

*/
  template<size_t N_CELL_FACETS>
    class MeshEntityArrays
    {
    
    private:
      dolfin::Mesh& _mesh;
      boost::multi_array<int, 2> cell_neighbors_array; // See docs above. This gets
                                                       // resized when needed.
      bool cell_neighbors_array_initialized{false};
      std::vector<std::array<double, N_CELL_FACETS*3>> cell_facet_normals_array; // See docs above.
      bool cell_facet_normals_array_initialized{false};
      
      static int NO_CELL; // Initialized in MeshEntityArrays.cpp
      
    public:

      //! The one and only explicit ctor.
      /*!

        \param mesh is the dolfin::Mesh for which the maps are computed.
        \param compute_field_mesh_maps is the flag to compute entity arrays that are used on the field mesh.
        \param compute_particle_mesh_maps is the flag to compute entity arrays that are used on the particle mesh.

      */
      MeshEntityArrays(dolfin::Mesh& mesh, bool compute_field_mesh_maps = false, bool compute_particle_mesh_maps = false);

      // The dtor (See pybind 8.5 Non-public destructors)
      ~MeshEntityArrays();

      // Compute an array containing the indices of neighboring cells for all cells
      void compute_cell_neighbors_array();
      // Obtain the cell neighors of one cell
      std::array<int, N_CELL_FACETS> get_cell_neighbors(std::size_t cell_index);

      // Compute cell facet-normals for all cells
      void compute_cell_facet_normals_array();
      // Obtain the facet-normals for one cell
      //      std::array<std::array<double, 3>, N_CELL_FACETS> get_cell_facet_normals(std::size_t cell_index);
      std::array<double, N_CELL_FACETS*3> get_cell_facet_normals(std::size_t cell_index);
      
      // Compute cell volumes
      // Probably not used on particle mesh?
    
      //    void compute_cell_volume_array(); // Compute cell volumes indexed by cell index

    };
  
} // namespace dnt

#endif
