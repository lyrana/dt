/*! \file particle_solib.cpp

  \brief This file creates a shared library containing the Python bindings for the
  C++ particle-movers in particle.h.

  This file contains the Python-to-C++ bindings allowing access to the C++
  particle-movers in particle.h.

  \namespace dnt

  \sa particle.h

*/
#include "particle.h"

// To create a .so file for a particular particle type, set the macro
// PARTICLE_MODULE_PREFIX in Makefile.part.

namespace py = pybind11;

// The PYBIND11_MODULE() macro creates a function that will be called when an import
// statement is issued from within Python. The module name is given as the first
// macro argument (it should not be in quotes). The second argument (m) defines a
// variable of type py::module which is the main interface for creating bindings. The
// method module::def() generates binding code that exposes the C++ functions to
// Python.

namespace dnt {

  PYBIND11_MODULE(PARTICLE_MODULE_NAME, m) {

    // C++ functions defined in particle.cpp
  
    // Interface to the C++ particle-advance functions

    m.def("move_charged_species_in_uniform_fields", &move_charged_species_in_uniform_fields<Ptype::PARTICLE_TYPE>);

    m.def("move_neutral_species_2_facets", &move_neutral_species<Ptype::PARTICLE_TYPE, 2>);
    m.def("move_neutral_species_3_facets", &move_neutral_species<Ptype::PARTICLE_TYPE, 3>);
    m.def("move_neutral_species_4_facets", &move_neutral_species<Ptype::PARTICLE_TYPE, 4>);
    
  
  } // ENDDEF: PYBIND11_MODULE()

} // namespace dnt
