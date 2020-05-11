/*! \file mesh_entity_arrays_solib.cpp

  \brief This file creates a shared library with the Python bindings for C++
         MeshEntityArrays objects, specialized for the number of facets in a mesh cell.

  This file contains the Python-to-C++ bindings allowing access to the C++ version of
  MeshEntityArrays from Python.

  \namespace dnt

  \sa MeshEntityArrays.h, MeshEntityArrays.cpp

*/
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include "MeshEntityArrays.cpp"

namespace py = pybind11;

namespace dnt {

  // The anonymous namespace limits the scope of the functions in it to this file.
  namespace {

    //! Make MeshEntityArrays classes for meshes with different numbers of cell facets.
    /*!

      makeMeshEntityArrays() is a 'helper' function. It causes the compiler to make a
      specialized MeshEntityArrays class for a mesh with the number of cell-facets
      specified by the N_CELL_FACETS template parameter. The class can then be
      instanced in C++ code. It can also be created from Python if the binding code
      is written below.

      \param N_CELL_FACETS a template parameter specifying the number of facets that
      a mesh-cell has.
      \param m is a py::module object created by PYBIND11_MODULE
      \param strFacets is a string used to create a unique Python class name.

      \return void

      \sa MeshEntityArrays.h, MeshEntityArrays.cpp

     */
    template <size_t N_CELL_FACETS>
    void makeMeshEntityArrays(py::module &m, std::string const & strFacets)
    {
      using MEA = MeshEntityArrays<N_CELL_FACETS>;
      // Make a class name with the particle structure type appended to
      // "MeshEntityArrays"
      std::string pyclass_name = std::string("MeshEntityArrays_") + strFacets;

      //      std::cout << "{mesh_entity_arrays_solib.cpp}makeMeshEntityArrays: pyclass_name " << pyclass_name << std::endl;
                                                      
      // Create the Python binding for a class with name pyclass_name.
      py::class_<MEA>(m, pyclass_name.c_str())

        // The ctor
        
        // The C++ ctor args types are template parameters in py::init<>()
        // In this case, the args are the dolfin mesh, and two boolean flags.
        // Note that if a py::arg() specification is needed (e.g., to specify a
        // default value), then every argument has to have a py::arg().
        .def(py::init<dolfin::Mesh&, bool, bool>(), py::arg("mesh"), py::arg("compute_field_mesh_maps") = false, py::arg("compute_particle_mesh_maps") = false)
        // The following statements create the bindings to MeshEntityArrays member
        // functions for a mesh with N_CELL_FACETS. The source code for these is in
        // MeshEntityArrays.cpp.  These aren't required to be callable from Python,
        // but can be used for convenient inspection of values(?)
        .def("compute_cell_facet_normals_array", &MEA::compute_cell_facet_normals_array)
        .def("compute_cell_neighbors_array", &MEA::compute_cell_neighbors_array)
        .def("get_cell_neighbors", &MEA::get_cell_neighbors)
        .def("get_cell_facet_normals", &MEA::get_cell_facet_normals);
      
      // cf. ~/workspace/dolfin/python/src/geometry.cpp for __getitem__, if needed.
        
    } // void makeMeshEntityArrays(py::module &m, std::string const & strFacets)

  } // namespace is anonymous
  

  // The PYBIND11_MODULE() macro creates a function that will be called when an import
  // statement is issued from within Python. The module name (mesh_entity_arrays_solib) is
  // given as the first macro argument (it should not be in quotes). The second argument (m)
  // defines a variable of type py::module which is the main interface for creating
  // bindings. The method module::def() generates binding code that exposes the C++
  // functions to Python.
  
  // Create a variable 'm' of type py::module
  PYBIND11_MODULE(mesh_entity_arrays_solib, m)
  {
    
    // Compile the classes "MeshEntityArrays_2/3/4_facets" and the Python bindings for
    // them.
    makeMeshEntityArrays<2>(m, "2_facets");
    makeMeshEntityArrays<3>(m, "3_facets");
    makeMeshEntityArrays<4>(m, "4_facets");
    
  } // ENDDEF: PYBIND11_MODULE(mesh_entity_arrays_solib, m)

  
} // namespace dnt
