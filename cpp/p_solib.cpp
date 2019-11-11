//#include "dolfin.h"
//#include "pseg.h"
#include "particle.h"

// To create a .so file for a particular type, set the type here, and in
// PARTICLE_MODULE_PREFIX. Then run 
// #define PARTICLE_TYPE cartesian_x
//#define PARTICLE_MODULE_PREFIX particle_pyb_
//#define PARTICLE_MODULE_NAME PARTICLE_MODULE_PREFIX ## PARTICLE_TYPE
// #define PARTICLE_MODULE_NAME p_pyb_cartesian_x

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
    //    m.def("move_neutral_species", &move_neutral_species<Ptype::PARTICLE_TYPE>);
    m.def("move_neutral_species", &move_neutral_species<Ptype::PARTICLE_TYPE, CELL_FACET_NUMBER>);
    
  
  } // ENDDEF: PYBIND11_MODULE()

} // namespace dnt
