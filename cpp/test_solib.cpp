/*! \file test_solib.cpp

  \brief This has the Python/C++ interface code for test.h

  This allows the C++ functions defined in test.h to be called from Python.
  Compile this using Makefile.test. The MODULE_NAME is given in the makefile.

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
// statement is issued from within Python. The module name (MODULE_NAME) is given as the
// first macro argument (it should not be in quotes). The second argument (m) defines
// a variable of type py::module which is the main interface for creating
// bindings. The method module::def() generates binding code that exposes the C++
// functions to Python.


namespace dnt {

  // Interface to C++ functions and types

  // Create a variable 'm' of type py::module.
  // MODULE_NAME can be specified using -DMODULE_NAME= in the makefile.  
  PYBIND11_MODULE(MODULE_NAME, m)
  {

    // C++ classes and functions declared and defined in test.h
    
    m.def("function_with_several_args", &function_with_several_args);
    m.def("function_with_DTcontrol_C_arg", &function_with_DTcontrol_C_arg);
    m.def("function_with_Particle_C_arg", &function_with_Particle_C_arg);
    m.def("function_with_pseg_arg", &function_with_pseg_arg<Ptype::PARTICLE_TYPE>);

  } // PYBIND11_MODULE(test_solib, m)

  
} // namespace dnt
