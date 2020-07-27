/*! \file numpy_types_solib.cpp

  \brief This file creates a shared library that registers types with Numpy

  In order to make Numpy arrays where the members are C structs, the struct and the
  variables it contains must be registered with Numpy.

  \namespace dnt

*/
// uncomment to disable assert()
// #define NDEBUG
#include <cassert>
#include "Pstruct.h"

namespace py = pybind11;

// The PYBIND11_MODULE() macro creates a function that will be called when an import
// statement is issued from within Python. The module name (dnt_cpp) is given as the
// first macro argument (it should not be in quotes). The second argument (m) defines
// a variable of type py::module which is the main interface for creating
// bindings. The method module::def() generates binding code that exposes the C++
// functions to Python.

namespace dnt {

  // Create a variable 'm' of type py::module.
  // MODULE_NAME can be specified using -DMODULE_NAME= in the makefile.  
  PYBIND11_MODULE(numpy_types_solib, m)
  {
    // The following registrations are needed for the Pstruct<> structured types to work with py::array_t. See pybind11: 12.2.3 Structured types.
    // The "_EX" variation allows the Python names of the variables in the structure to be different from the variable names in the C++ struct.
    PYBIND11_NUMPY_DTYPE_EX(Pstruct<Ptype::cartesian_x>, x_, "x", x0_, "x0", ux_, "ux", weight_, "weight", bitflags_, "bitflags", cell_index_, "cell_index", unique_ID_, "unique_ID", crossings_, "crossings");
    PYBIND11_NUMPY_DTYPE_EX(Pstruct<Ptype::cartesian_xy>, x_, "x", y_, "y", x0_, "x0", y0_, "y0", ux_, "ux", uy_, "uy", weight_, "weight", bitflags_, "bitflags", cell_index_, "cell_index", unique_ID_, "unique_ID", crossings_, "crossings");
    PYBIND11_NUMPY_DTYPE_EX(Pstruct<Ptype::cartesian_xyz>, x_, "x", y_, "y", z_, "z", x0_, "x0", y0_, "y0", z0_, "z0", ux_, "ux", uy_, "uy", uz_, "uz", weight_, "weight", bitflags_, "bitflags", cell_index_, "cell_index", unique_ID_, "unique_ID", crossings_, "crossings");

  } // ENDDEF: PYBIND11_MODULE(numpy_types_solib, m)
  
} // namespace dnt
