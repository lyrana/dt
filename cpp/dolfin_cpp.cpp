/*! \file dolfin_cpp.cpp

  \brief This file creates a shared library with the Python bindings for
  the C++ dolfin-based functions.

  \namespace dnt
  \sa dolfin.h, dolfin.cpp, segmentedarraypair_cpp.cpp

  This file contains the Python-to-C++ bindings allows access to the C++ versions of functions
  defined in dolfin.cpp.

*/

#include "dolfin.h"

namespace py = pybind11;

// The PYBIND11_MODULE() macro creates a function that will be called when an import
// statement is issued from within Python. The module name (dnt_cpp) is given as the
// first macro argument (it should not be in quotes). The second argument (m) defines
// a variable of type py::module which is the main interface for creating
// bindings. The method module::def() generates binding code that exposes the C++
// functions to Python.


namespace dnt {

  // Helper functions for TEMPLATED CLASSES are put into an anonymous namespace,
  // which makes them available only in this file. They're used to create
  // specializations of the classes (see bottom of this file).
  
  namespace {
    /*    
    // Interface to the C++ class SegmentedArrayPair_Cpp
    template <typename PS>
    void makeSegmentedArrayPair_Cpp(py::module &m, std::string const & PStype)
    {
    using SAP = SegmentedArrayPair_Cpp<PS>;
    std::string pyclass_name = std::string("SegmentedArrayPair_Cpp") + PStype;
    py::class_<SAP>(m, pyclass_name.c_str())

    // ctors
    // The ctor args types become template parameters
    //        .def(py::init<int>())
    .def(py::init<py::ssize_t>())

    // member functions
    .def("get_number_of_segments", &SAP::get_number_of_segments)
    .def("put", &SAP::put);
    }
    */

    /*
      void declare_array(py::module &m, std::string &typestr) {
      using Class = Array2D<T>;
      std::string pyclass_name = std::string("Array2D") + typestr;
      py::class_<Class>(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr())
      .def(py::init<>())
      .def(py::init<Class::xy_t, Class::xy_t, T>())
      .def("size",      &Class::size)
      .def("width",     &Class::width)
      .def("height",    &Class::height);
    */
      
  } // anonymous namespace

  // Create a variable 'm' of type py::module
  PYBIND11_MODULE(dolfin_cpp, m) {

    m.def("cell_contains_point", &cell_contains_point);

    // Call the C++ functions defined in dolfin.cpp
  
    //  void dolfin_DnT_pstruct_instances()
    dolfin::Mesh mesh;
    dolfin::Mesh& meshRef = mesh;
    py::array_t<int> vertices;
    double* point(nullptr);

    // Do these calls need to be inside the MODULE?
    cell_contains_point(meshRef, vertices, point);
  }

} // namespace dnt
