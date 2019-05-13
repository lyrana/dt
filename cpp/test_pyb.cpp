/*! \file test_pyb.cpp

  \brief This has the Python/C++ interface code for test.h

  \namespace dnt

*/
#include "test.h"

namespace py = pybind11;

// The PYBIND11_MODULE() macro creates a function that will be called when an import
// statement is issued from within Python. The module name (test_cpp) is given as the
// first macro argument (it should not be in quotes). The second argument (m) defines
// a variable of type py::module which is the main interface for creating
// bindings. The method module::def() generates binding code that exposes the C++
// functions to Python.


namespace dnt {

  // Interface to C++ functions and types

  // Create a variable 'm' of type py::module
  PYBIND11_MODULE(test_pyb, m)
  {

// Register Pstruct<cartesian_x> as a Numpy dtype descriptor. It is of type py::dtype (a pybind11 class).
// The "_EX" variation allows the Python names of the variables in the structure to be different from the variable names in the C++ struct.    
    PYBIND11_NUMPY_DTYPE_EX(Pstruct<Ptype::cartesian_x>, x_, "x", x0_, "x0", ux_, "ux", weight_, "weight", bitflags_, "bitflags", cell_index_, "cell_index", unique_ID_, "unique_ID", crossings_, "crossings");

// Register 2D and 3D cartesian types
    PYBIND11_NUMPY_DTYPE_EX(Pstruct<Ptype::cartesian_x_y>, x_, "x", y_, "y", x0_, "x0", y0_, "y0", ux_, "ux", uy_, "uy", weight_, "weight", bitflags_, "bitflags", cell_index_, "cell_index", unique_ID_, "unique_ID", crossings_, "crossings");
    PYBIND11_NUMPY_DTYPE_EX(Pstruct<Ptype::cartesian_x_y_z>, x_, "x", y_, "y", z_, "z", x0_, "x0", y0_, "y0", z0_, "z0", ux_, "ux", uy_, "uy", uz_, "uz", weight_, "weight", bitflags_, "bitflags", cell_index_, "cell_index", unique_ID_, "unique_ID", crossings_, "crossings");

    // C++ classes and functions declared and ?defined in test.h

    m.def("function_with_DTcontrol_arg", &function_with_DTcontrol_arg);
    
    // Templated functions need to be created with the particular types needed specified
    m.def("function_with_SAP_cartesian_x_arg", &function_with_SAP_arg<Ptype::cartesian_x>);
    m.def("function_with_SAP_cartesian_x_y_arg", &function_with_SAP_arg<Ptype::cartesian_x_y>);
    m.def("function_with_SAP_cartesian_x_y_z_arg", &function_with_SAP_arg<Ptype::cartesian_x_y_z>);
    
  } // PYBIND11_MODULE(test_pyb, m)

  
} // namespace dnt
