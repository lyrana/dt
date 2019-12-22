/*! \file MeshEntityArrays.cpp

  \brief This file defines the C++ class MeshEntityArrays.

    From MeshEntity.h, the definition of entities() is:
      Return array of indices for incident mesh entities of given
      topological dimension

      @param     dim (std::size_t)
                 The topological dimension.

      @return     std::size_t
                  The index for incident mesh entities of given dimension.
      const unsigned int* entities(std::size_t dim) const
      {
        const unsigned int* initialized_mesh_entities = _mesh->topology()(_dim, dim)(_local_index);
        dolfin_assert(initialized_mesh_entities);
        return initialized_mesh_entities; // The returned array of indices.
    }

  Connectivity is computed in TopologyComputation::compute_connectivity()


  \namespace dnt
  \sa MeshEntityArrays.h

*/

#include "MeshEntityArrays.h"
//#include "predicates.h"

//namespace py = pybind11;

/*
void init_ex2(py::module &m) {
m.def("sub", [](int a, int b) { return a - b; });
}
*/

namespace dnt
{

  template<size_t N_CELL_FACETS>
  int MeshEntityArrays<N_CELL_FACETS>::NO_CELL = -1;
    
  // BEGINDEF: MeshEntityArrays ctor
  //! The one and only MeshEntityArrays ctor.
  /*!

    \param mesh is the dolfin::Mesh for which the maps are computed.
    \param compute_field_mesh_maps is the flag to compute entity arrays that are used on the field mesh.
    \param compute_particle_mesh_maps is the flag to compute entity arrays that are used on the particle mesh.

  */
  template<size_t N_CELL_FACETS>
  MeshEntityArrays<N_CELL_FACETS>::MeshEntityArrays(dolfin::Mesh& mesh, bool compute_field_mesh_maps, bool compute_particle_mesh_maps):
    _mesh(mesh)
  {
    //    std::cout << "Hello from the MeshEntityArrays constructor" << std::endl;

    if (compute_field_mesh_maps == true)
      {
        //
      }
    
    if (compute_particle_mesh_maps == true)
      {
        //        std::cout << "Computing cell_facet_normals_array" << std::endl;
        compute_cell_facet_normals_array();
        //        std::cout << "Computing cell_neighbors_array" << std::endl;
        compute_cell_neighbors_array();
      }
  }
  // ENDDEF: MeshEntityArrays ctor

  // The dtor (See pybind 8.5 Non-public destructors)
  template<size_t N_CELL_FACETS>
  MeshEntityArrays<N_CELL_FACETS>::~MeshEntityArrays()
  {
    std::cout << "The dtor ~MeshEntityArrays has been called" << std::endl;
      
    // Release the arrays?
      
  } 
  
  // BEGINDEF: compute_cell_neighbors_array
  //! Make a 2D array that gives the indices of cells that share a common facet.
  /*!

    The 2D array cell_neighbors_array[i][j] is indexed by (1) cell index i (2)
    local facet number j.

    cell_neighbors_array[i] is a list of the indices of cells attached to cell i.

    cell_neighbors_array[i][j] is the index of the cell attached to facet j of cell i.
        
  */
  template<size_t N_CELL_FACETS>
  void MeshEntityArrays<N_CELL_FACETS>::compute_cell_neighbors_array()
  {

    //    std::map <int, int*> cell_neighbors_array;
    
    // Use "tDim" to compute topological quantities (e.g., number of facets in a cell).
    auto tDim = _mesh.topology().dim();

    // Generate facet-to-cell connectivity data, i.e., given a mesh-level facet index,
    // list the cells on either side.
    _mesh.init(tDim - 1, tDim);

    // For every cell, build a list of cells that are connected to its facets
    // but are not the iterated cell.

    // Double loop: outer one over cells in the mesh, inner one over facets of the cell.
    
    //    auto nFacets = tDim + 1;

    //    std::cout << "Resizing the cell_neighbors_array for " << _mesh.num_cells() << " cells" << std::endl;

    cell_neighbors_array.resize(boost::extents[_mesh.num_cells()][N_CELL_FACETS]);
    
    // Create the list of cells sharing a facet with current cell    
    for (dolfin::CellIterator cIter(_mesh); !cIter.end(); ++cIter)
      {
        // Make an array to hold the indices of the neighbor cells.  A cell has
        // nFacets, and each facet can have at most 1 neighbor-cell attached.
        
        //        auto neighborCells = new int[nFacet]; // Room for the indices of neighbor cells. This has to be persistent storage.
        //        auto thisCellIndex = *cIter.index();
        auto thisCellIndex = cIter->index();

        //        std::cout << "thisCellIndex " << thisCellIndex << std::endl;

        auto neighborCells = cell_neighbors_array[thisCellIndex];

        //        std::cout << "Entering the facet loop" << std::endl;
        
        // Loop on the facets of this cell
        //        for (dolfin::FacetIterator fIter(*cIter), int i_facet = 0; !fIter.end(); ++fIter, ++i_facet)
        size_t fi = 0;
        for (dolfin::FacetIterator fIter(*cIter); !fIter.end(); ++fIter)
          {
            // Obtain the number of cells attached to this facet. There are at most 2: the
            // current cell and the cell on the other side (unless the facet is on the
            // boundary).
            auto numberOfAttachedCells = fIter->num_entities(tDim);
            //            std::cout << "numberOfAttachedCells " << numberOfAttachedCells << std::endl;
            assert(numberOfAttachedCells >= 1 && numberOfAttachedCells <= 2);

            auto attachedCellIndices = fIter->entities(tDim);
            /*            
            for (size_t ic = 0; ic < numberOfAttachedCells; ic++)
              {
                if (attachedCellIndices[ic] != thisCellIndex)
                  {
                    std::cout << "attachedCellIndices[ic] " << attachedCellIndices[ic] << std::endl;
                  }
              }
            */
            
            if (numberOfAttachedCells == 2)
              {
                // If the first one is the current cell, the neighbor is the second one.
                if (attachedCellIndices[0] == thisCellIndex)
                  {
                    neighborCells[fi] = attachedCellIndices[1];
                  }
                else // If not, the neighbor is the first one
                  {
                    neighborCells[fi] = attachedCellIndices[0];
                  }
              }
            else // There's only one cell attached to the facet. It must be on a boundary.
              {
                neighborCells[fi] = NO_CELL;
              }
            fi++;
          }
      } // for (dolfin::CellIterator cIter(_mesh); !cIter.end(); ++cIter)

    // Check the neighbors of cell 0 (Works only for 3-facet cells, i.e., triangles)
    
    // std::cout << "The cell indices of the neighbors of cell 0 are " << get_cell_neighbors(0)[0] << " " << get_cell_neighbors(0)[1] << " " << get_cell_neighbors(0)[2] << std::endl;

    cell_neighbors_array_initialized = true;
    
    return;
  }
  // ENDDEF: compute_cell_neighbors_array


  // BEGINDEF: get_cell_neighbors
  //! Return the list of cell indices of cells that share a facet with the current cell.
  /*!  

    \return The neighbor-cell indices in std::array<int, N_CELL_FACETS> neighborCells

  */
  template<size_t N_CELL_FACETS>
  std::array<int, N_CELL_FACETS> MeshEntityArrays<N_CELL_FACETS>::get_cell_neighbors(std::size_t cell_index)
  {
    // int* neighborCellsArrayPtr = cell_neighbors_array.data();

    std::array<int, N_CELL_FACETS> neighborCells;
    
    if (cell_neighbors_array_initialized)
      {
        for (size_t fi = 0; fi < N_CELL_FACETS; fi++)
          {
            neighborCells[fi] = cell_neighbors_array[cell_index][fi];
          }
      }
    else
      {
        std::cout << "{MeshEntityArrays.cpp}get_cell_neighbors(): cell_neighbors_array has not been initialized" << std::endl;
        std::exit(EXIT_FAILURE);        
      }

    return neighborCells;
  }
  // ENDDEF: get_cell_neighbors

  
  // BEGINDEF: compute_cell_facet_normals_array
  //! Make a 1D array giving the cell facet-normal 3-vectors, indexed by the cell index.
  /*!  
    The normals are unit vectors represented by Points (3-vectors; see below).  For
    tetrahedra, they are computed in dolfin/dolfin/mesh/TetrahedronCell.cpp.

    From dolfin/dolfin/mesh/Cell.h:
    Compute normal of given facet with respect to the cell
    Point normal(std::size_t facet) const
    { return _mesh->type().normal(*this, facet); }

    From dolfin/dolfin/geometry/Point.h:
    A Point represents a point in :math:`\mathbb{R}^3` with
    coordinates :math:`x, y, z,` or alternatively, a vector in
    :math:`\mathbb{R}^3`, supporting standard operations like the
    norm, distances, scalar and vector products etc.

    cell_facet_normals has type: std::vector<std::array<double, N_CELL_FACETS*3>>
    i.e., it's a std:vector with length equal to to number of cells in the mesh. Each
    item in the array is an std::array of doubles. For example:
        cell_facet_normals_array[i] contains the array [n0, n1, ...], where n0, n1,
        ... are the unit normal 3-vectors on each facet of cell i.
    The normal vector always has 3 components, regardless of the dimension of the mesh.

    \return void

  */
  template<size_t N_CELL_FACETS>
  void MeshEntityArrays<N_CELL_FACETS>::compute_cell_facet_normals_array()
  {
    
    //    const std::size_t tDim = _mesh.topology().dim();
    //    const std::size_t gDim = _mesh.geometry().dim();
    //    size_t nFacets = tDim + 1;

    cell_facet_normals_array.resize(_mesh.num_cells());
    
    for (dolfin::CellIterator cIter(_mesh); !cIter.end(); ++cIter)
      {
        auto thisCellIndex = cIter->index();
        //        std::cout << "thisCellIndex " << thisCellIndex << std::endl;        
        std::array<double, N_CELL_FACETS*3> normals;

        //std::cout << "normals type is " << typeid(normals).name() << std::endl;
   
        for (size_t fi = 0; fi < N_CELL_FACETS; fi++)
          {
            normals[3*fi]   = cIter->normal(fi).x();
            normals[3*fi+1] = cIter->normal(fi).y();
            normals[3*fi+2] = cIter->normal(fi).z();
          }
        cell_facet_normals_array[thisCellIndex] = normals;
      }

    //    std::cout << "cell_facet_normals_array[0] is " << cell_facet_normals_array[0][0] << " " << cell_facet_normals_array[0][1] << " " << cell_facet_normals_array[0][2] << std::endl;
    
    cell_facet_normals_array_initialized = true;
    
    return;
  }
  // ENDDEF: void compute_cell_facet_normals_array


  // BEGINDEF: get_cell_facet_normals
  //! Return the list of cell indices of cells that share a facet with the current cell.
  /*!  

    \return The facet-normal vectors in a std::array<std::array<double, 3>, N_CELL_FACETS>

  */
  template<size_t N_CELL_FACETS>
    //  std::array<std::array<double, 3>, N_CELL_FACETS> MeshEntityArrays<N_CELL_FACETS>::get_cell_facet_normals(std::size_t cell_index)
  //    std::array<double, N_CELL_FACETS*3> MeshEntityArrays<N_CELL_FACETS>::get_cell_facet_normals(std::size_t cell_index)
    std::array<double, N_CELL_FACETS*3>& MeshEntityArrays<N_CELL_FACETS>::get_cell_facet_normals(std::size_t cell_index)
    
  {
    //    std::cout << "Entered get_cell_facet_normals()" << std::endl;
    
    //    std::array<double, N_CELL_FACETS*3> cellFacetNormals;
    //    std::array<std::array<double, 3>, N_CELL_FACETS> cellFacetNormals; // This messes up the storage order

    //    auto computedNormals1 = cell_facet_normals_array[cell_index];
    //    std::cout << "get_cell_facet_normals() [0] is " << computedNormals1[0] << std::endl;
    
    if (cell_facet_normals_array_initialized)
      {
        //        auto cellFacetNormals = cell_facet_normals_array[cell_index];
        //        auto computedNormals = cell_facet_normals_array[cell_index];
        /*
        for (size_t fi = 0; fi < N_CELL_FACETS; fi++)
          {
            cellFacetNormals[0][fi] = computedNormals[3*fi];
            cellFacetNormals[1][fi] = computedNormals[3*fi+1];
            cellFacetNormals[2][fi] = computedNormals[3*fi+2];
          }
        */
      }
    else
      {
        std::cout << "{MeshEntityArrays.cpp}get_cell_facet_normals(): cell_facet_normals_array has not been initialized" << std::endl;
        std::exit(EXIT_FAILURE);
      }

    //    return cellFacetNormals;
    return cell_facet_normals_array[cell_index];
  }
  // ENDDEF: get_cell_facet_normals

} // namespace dnt

template class dnt::MeshEntityArrays<2>;
template class dnt::MeshEntityArrays<3>;
template class dnt::MeshEntityArrays<4>;


