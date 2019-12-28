/*! \file test_solib.cpp

  \brief This has the Python/C++ interface code for test.h

  \namespace dnt

*/

#include <memory>
#include <iostream>

#include "SegmentedArrayPair.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
// #include <pybind11/eigen.h>
#include <pybind11/operators.h>

#include "test.h"

namespace py = pybind11;

// The PYBIND11_MODULE() macro creates a function that will be called when an import
// statement is issued from within Python. The module name (test_solib) is given as the
// first macro argument (it should not be in quotes). The second argument (m) defines
// a variable of type py::module which is the main interface for creating
// bindings. The method module::def() generates binding code that exposes the C++
// functions to Python.


namespace dnt {

  // Interface to C++ functions and types

  // Create a variable 'm' of type py::module
  PYBIND11_MODULE(MODULE_NAME, m)
  {

    // C++ classes and functions declared and defined in test.h
    m.def("function_with_DTcontrol_arg", &function_with_DTcontrol_arg);
    m.def("function_with_particle_P_arg", &function_with_particle_P_arg);

  } // PYBIND11_MODULE(test_solib, m)

  
} // namespace dnt
