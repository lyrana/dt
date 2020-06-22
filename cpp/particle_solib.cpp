/*! \file particle_solib.cpp

  \brief This file creates a shared library containing the Python bindings for the
  C++ particle-movers in particle.h.

  This file contains the Python-to-C++ bindings allowing access to the C++
  particle-movers in particle.h.

  \namespace dnt

  \sa particle.h

*/
#include "particle.h"

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
// Using intermediate functions causes concatenation to be delayed, so that names are
// expanded before concatenation
#define ACSIEF(x) advance_charged_species_in_E_field_##x
#define ADVANCE_CHARGED_SPECIES_IN_E_FIELD_(x) ACSIEF(x)
#define ANS(x) advance_neutral_species_##x
#define ADVANCE_NEUTRAL_SPECIES_(x) ANS(x)

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

  // Set the bit patterns for flags. These are static class members, so they're set
  // outside the ctor, in a C++ source file, so they're compiled only once. (If they're in
  // a header file that's included in more than one file, that generates an error)
  //int Pstruct<Ptype::PARTICLE_TYPE>::DELETE_FLAG = 0b1;  // the lowest bit is 1
  //int Pstruct<Ptype::PARTICLE_TYPE>::TRAJECTORY_FLAG = 0b1 << 1; // the second lowest bit  // Need all types because the source-code in particle.h has them all
  int Pstruct<Ptype::cartesian_xyz>::DELETE_FLAG = 0b1;  // the lowest bit is 1
  int Pstruct<Ptype::cartesian_xyz>::TRAJECTORY_FLAG = 0b1 << 1; // the second lowest bit 
  int Pstruct<Ptype::cartesian_xy>::DELETE_FLAG = 0b1;  // the lowest bit is 1
  int Pstruct<Ptype::cartesian_xy>::TRAJECTORY_FLAG = 0b1 << 1; // the second lowest bit 
  int Pstruct<Ptype::cartesian_x>::DELETE_FLAG = 0b1;  // the lowest bit is 1
  int Pstruct<Ptype::cartesian_x>::TRAJECTORY_FLAG = 0b1 << 1; // the second lowest bit 
  
  PYBIND11_MODULE(MODULE_NAME, m) {

    // Interface to the C++ particle-advance functions
    // Note that PARTICLE_TYPE here is defined in Makefile.part

    // From https://en.cppreference.com/w/cpp/language/function_template:
    // A function template by itself is not a type, or a function, or any other entity. No
    // code is generated from a source file that contains only template definitions. In order
    // for any code to appear, a template must be instantiated: the template arguments must be
    // determined so that the compiler can generate an actual function (or class, from a class
    // template).
    
    m.def("initialize_particle_integration", &initialize_particle_integration);    

    // Interface to the ParticleMeshBoundaryConditions class is defined in particle.h.
    makeParticleMeshBoundaryConditions<Ptype::PARTICLE_TYPE>(m,  TOSTRING(PARTICLE_TYPE));

    // xyz version of .so

    // m.def("advance_charged_species_in_uniform_fields", &advance_charged_species_in_uniform_fields_cartesian_xyz);
    // m.def("advance_neutral_species_4_facets", &advance_neutral_species_cartesian_xyz<4>);
    // m.def("advance_neutral_species_3_facets", &advance_neutral_species_cartesian_xyz<3>);
    // m.def("advance_neutral_species_2_facets", &advance_neutral_species_cartesian_xyz<2>);
    

    // These calls create needed specializations

    //makeParticleMeshBoundaryConditions<Ptype::cartesian_xyz>(m, "cartesian_xyz");

    
    // xy version of .so
    //m.def("advance_charged_species_in_E_field_3_facets", &advance_charged_species_in_E_field_cartesian_xy<3>, py::arg("particle_P"), py::arg("species_name"), py::arg("ctrl"), py::arg("neg_E_field") = nullptr, py::arg("external_E_field") = nullptr, py::arg("accel_only") = false);
    m.def("advance_charged_species_in_E_field_3_facets", &ADVANCE_CHARGED_SPECIES_IN_E_FIELD_(PARTICLE_TYPE)<3>, py::arg("particle_P"), py::arg("species_name"), py::arg("ctrl"), py::arg("neg_E_field") = nullptr, py::arg("external_E_field") = nullptr, py::arg("accel_only") = false);
    // m.def("advance_neutral_species_3_facets", &advance_neutral_species_cartesian_xy<3>);
    m.def("advance_neutral_species_3_facets", &ADVANCE_NEUTRAL_SPECIES_(PARTICLE_TYPE)<3>);
    // Interface to the ParticleMeshBoundaryConditions class is defined in particle.h.
    // These calls create needed specializations
    // makeParticleMeshBoundaryConditions<Ptype::cartesian_xy>(m, "cartesian_xy");

  } // ENDDEF: PYBIND11_MODULE(MODULE_NAME, m)


} // namespace dnt
