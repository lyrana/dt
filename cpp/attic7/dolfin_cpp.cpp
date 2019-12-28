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




// // Create a variable 'm' of type py::module
PYBIND11_MODULE(dolfin_cpp, m) {

// typeinfo may be registered before the dtype descriptor for scalar casts to work...

///
// Allow Python types to be used in C++
///  
  
//  py::class_() creates Python bindings for a C++ class or struct-style data
// structure.  The following statements connect the Python names DnT_pstruct?D to C++
// structs with the same names.
  py::class_<DnT_pstruct1D>(m, "DnT_pstruct1D");
  py::class_<DnT_pstruct2D>(m, "DnT_pstruct2D");
  py::class_<DnT_pstruct3D>(m, "DnT_pstruct3D");

// The following statements allow the Numpy structured types DnT_pstructXD to be used in C++.
// In particular they can be used as template arguments to the Numpy array type py::array_t.
  
// Register DnT_pstruct1D as a Numpy dtype descriptor. DnT_pstruct1D is of type py::dtype (a pybind11 class).
//  PYBIND11_NUMPY_DTYPE(DnT_pstruct1D, x, x0, ux, weight, bitflags, cell_index, unique_ID, crossings);
// The "_EX" variation allows the Python names of the variables in the structure to be different from the variable names in the C++ struct.
  PYBIND11_NUMPY_DTYPE_EX(DnT_pstruct1D, x_, "x", x0_, "x0", ux_, "ux", weight_, "weight", bitflags_, "bitflags", cell_index_, "cell_index", unique_ID_, "unique_ID", crossings_, "crossings");

// Register DnT_pstruct2D and DnT_pstruct3D
  PYBIND11_NUMPY_DTYPE_EX(DnT_pstruct2D, x_, "x", y_, "y", x0_, "x0", y0_, "y0", ux_, "ux", uy_, "uy", weight_, "weight", bitflags_, "bitflags", cell_index_, "cell_index", unique_ID_, "unique_ID", crossings_, "crossings");
  PYBIND11_NUMPY_DTYPE_EX(DnT_pstruct3D, x_, "x", y_, "y", z_, "z", x0_, "x0", y0_, "y0", z0_, "z0", ux_, "ux", uy_, "uy", uz_, "uz", weight_, "weight", bitflags_, "bitflags", cell_index_, "cell_index", unique_ID_, "unique_ID", crossings_, "crossings");

  
///
// Bindings for mesh operations
///  

// dolfin.cpp:

// Connect the Python symbol divide_by_cell_volumes() to the C++ function declared as
//         void divide_by_cell_volumes(dolfin::Function& dF, py::dict dict)
  m.def("divide_by_cell_volumes", &divide_by_cell_volumes);
  
///
// Bindings for particle-to-mesh operations
///  

// pseg.cpp:
  
// Connect the Python symbol add_weights_to_cells() to the C++ function declared as
//       void add_weights_to_cells(py::array_t<S, 0> pseg, dolfin::Function& dF) {}
  
// These use a full cast of the function arguments:
//  m.def("add_weights_to_cells1D", (void (*) (pybind11::array_t<DnT_pstruct1D, 0>, dolfin::Function&)) &add_weights_to_cells);
//  m.def("add_weights_to_cells2D", (void (*) (pybind11::array_t<DnT_pstruct2D, 0>, dolfin::Function&)) &add_weights_to_cells);
//  m.def("add_weights_to_cells3D", (void (*) (pybind11::array_t<DnT_pstruct3D, 0>, dolfin::Function&)) &add_weights_to_cells);
  
// This simpler binding also makes the compiler generate the special cases:
  m.def("add_weights_to_cells1D", &add_weights_to_cells<DnT_pstruct1D>);
  m.def("add_weights_to_cells2D", &add_weights_to_cells<DnT_pstruct2D>);
  m.def("add_weights_to_cells3D", &add_weights_to_cells<DnT_pstruct3D>);

// Connect the Python symbol interpolate_weights_to_dofs to the C++ function declared as
//       void interpolate_weights_to_dofs(py::array_t<S, 0> pseg, dolfin::Function& dF) {}

// These use a full cast of the function arguments:  
//  m.def("interpolate_weights_to_dofs1D", (void (*) (pybind11::array_t<DnT_pstruct1D, 0>, dolfin::Function&)) &interpolate_weights_to_dofs);
//  m.def("interpolate_weights_to_dofs2D", (void (*) (pybind11::array_t<DnT_pstruct2D, 0>, dolfin::Function&)) &interpolate_weights_to_dofs);
//  m.def("interpolate_weights_to_dofs3D", (void (*) (pybind11::array_t<DnT_pstruct3D, 0>, dolfin::Function&)) &interpolate_weights_to_dofs);
  
// This simpler binding also makes the compiler generate the special cases:
  m.def("interpolate_weights_to_dofs1D", &interpolate_weights_to_dofs<DnT_pstruct1D>);
  m.def("interpolate_weights_to_dofs2D", &interpolate_weights_to_dofs<DnT_pstruct2D>);
  m.def("interpolate_weights_to_dofs3D", &interpolate_weights_to_dofs<DnT_pstruct3D>);
  
///
// Bindings for mesh-to-particle functions
///


///
// Bindings for printing functions
///


///
// Bindings for other functions
///

// C++ functions defined in dolfin.cpp:

// Works with std=c++11:  
  m.def("cell_contains_point", (bool (*) (dolfin::Mesh&, py::array_t<int>, Pstruct<Ptype::cartesian_x>)) &cell_contains_point);
  m.def("cell_contains_point", (bool (*) (dolfin::Mesh&, py::array_t<int>, Pstruct<Ptype::cartesian_xy>)) &cell_contains_point);
  m.def("cell_contains_point", (bool (*) (dolfin::Mesh&, py::array_t<int>, Pstruct<Ptype::cartesian_xyz>)) &cell_contains_point);

  m.def("cell_contains_point", (bool (*) (dolfin::Mesh&, py::array_t<int>, double&)) &cell_contains_point);
  m.def("cell_contains_point", (bool (*) (dolfin::Mesh&, py::array_t<int>, double&, double&)) &cell_contains_point);
  m.def("cell_contains_point", (bool (*) (dolfin::Mesh&, py::array_t<int>, double&, double&, double&)) &cell_contains_point);
  
// Works with std=c++14:  
//  m.def("cell_contains_point", py::overload_cast<dolfin::Mesh&, py::array_t<int>, DnT_pstruct1D>(&cell_contains_point), "1D point");
//  m.def("cell_contains_point", py::overload_cast<dolfin::Mesh&, py::array_t<int>, DnT_pstruct2D>(&cell_contains_point), "2D point");
//  m.def("cell_contains_point", py::overload_cast<dolfin::Mesh&, py::array_t<int>, DnT_pstruct3D>(&cell_contains_point), "3D point");

// Example: Connect the Python symbol f_simple() to the C++ function in braces {...}
// This prints only one record of the array
  m.def("f_simple", [](DnT_pstruct1D p) { return p.x_ * 10; });


// Call the C++ functions defined in dolfin.cpp
  

  //  void dolfin_DnT_pstruct_instances()
{
  dolfin::Mesh mesh;
  dolfin::Mesh& meshRef = mesh;
  py::array_t<int> vertices;
  Pstruct<Ptype::cartesian_x> ps_x;
  Pstruct<Ptype::cartesian_xy> ps_xy;
  Pstruct<Ptype::cartesian_xyz> ps_xyz;

  // c++filt said that the following symbols were not in the .so file.  These
  // statements force the compiler to make these specialized versions of the
  // templated function in dolfin.o.

  // Do these calls need to be inside the MODULE?
  cell_contains_point(meshRef, vertices, ps1D);
  cell_contains_point(meshRef, vertices, ps2D);
  cell_contains_point(meshRef, vertices, ps3D);
}

}

} // namespace dnt
