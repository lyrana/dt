/*! \file user_particle_boundary_functions_solib.cpp

  \brief This file creates a shared library containing the Python bindings for particle call-back functions in UserParticleBoundaryFunctions_2D_e.h.

  This file contains the Python-to-C++ bindings allowing Python to access the C++
  call-back functions in UserParticleBoundaryFunctions_2D_e.h. These functions are usually
  called directly from C++, but they are referenced in a Python dictionary in a
  ParticleMeshBoundaryConditions_C instance, which associates boundaries with
  call-back functions.

  \namespace dnt

  \sa particle.h ParticleMeshBoundaryConditions.h UserParticleBoundaryFunctions_2D_e.h

*/
#include "UserParticleBoundaryFunctions_2D_e.h"

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

namespace py = pybind11;

namespace dnt {

  // Set the bit patterns for flags. These are static class members, so they're set
  // outside the ctor, in a C++ source file, so they're compiled only once. (If they're in
  // a header file that's included in more than one file, that generates an error)
  // This variable is used in UserParticleBoundaryFunctions_2D_e.h:
  int Pstruct<Ptype::PARTICLE_TYPE>::DELETE_FLAG = 0b1;  // the lowest bit is 1
  
  // The anonymous namespace limits the scope of the functions in it to this file.
  namespace {

    //! Make UserParticleBoundaryFunctions classes to treat particles of different Ptypes.
    /*!

      makeUserParticleBoundaryFunctions() is a 'helper' function. It makes the Python
      bindings for a specialized UserParticleBoundaryFunctions class and its member
      functions based on the PT template parameter. Class instances can then be created
      and used in Python.

      \param PT is a template parameter for the function, specifying the class Ptype.
      \param m is a py::module object created by PYBIND11_MODULE.
      \param PT_str is a string used to create a unique Python class name based on the Ptype.

      \return void

      \sa UserParticleBoundaryFunctions_2D_e.h, UserParticleBoundaryFunctions_2D_e.cpp

     */
    template <Ptype PT>
    void makeUserParticleBoundaryFunctions(py::module &m, std::string const & PT_str)
    {
      using UPBF = UserParticleBoundaryFunctions<PT>;
      // Make a class name with the particle structure type (Ptype) string appended to the
      // string "UserParticleBoundaryFunctions_"
      
      // std::string pyclass_name = std::string("UserParticleBoundaryFunctions_") + PT_str;
      // Don't need the PT_str in the name since the MODULE_NAME has it:
      std::string pyclass_name = std::string("UserParticleBoundaryFunctions");
      
      // Create the Python binding for this class
      py::class_<UPBF>(m, pyclass_name.c_str())
        
        // The ctor is UserParticleBoundaryFunctions(position_coordinates, dx)
        //    :param position_coordinates: Example: ['x', 'y',]
        //    :vartype position_coordinates: list of strings
        // The C++ ctor args types are template parameters in py::init<>()
        // Note that if a py::arg() specification is needed (e.g., to specify a
        // default value), then every argument has to have a py::arg().
        .def(py::init< py::list& >())

        // The following statements create the bindings to UserParticleBoundaryFunctions member
        // functions for a mesh with Ptype PT.
        // Don't seem to need this.
        .def("default_bc", &UPBF::default_bc);
      
      // cf. ~/workspace/dolfin/python/src/geometry.cpp for __getitem__, if needed.
        
    } // void makeUserParticleBoundaryFunctions

  } // namespace is anonymous

      
// The PYBIND11_MODULE() macro creates a function that will be called when an import
// statement is issued from within Python. The module name is given as the first
// macro argument (it should not be in quotes). The second argument (m) defines a
// variable of type py::module which is the main interface for creating bindings. The
// method module::def() generates binding code that exposes the C++ functions to
// Python.

// Create a variable 'm' of type py::module  
// MODULE_NAME can be specified using -DMODULE_NAME= in the makefile.
  PYBIND11_MODULE(MODULE_NAME, m) {

    // Interface to the C++ particle call-back functions.
    // Note that PARTICLE_TYPE here is defined in Makefile.part

    // From https://en.cppreference.com/w/cpp/language/function_template:
    // A function template by itself is not a type, or a function, or any other entity. No
    // code is generated from a source file that contains only template definitions. In order
    // for any code to appear, a template must be instantiated: the template arguments must be
    // determined so that the compiler can generate an actual function (or class, from a class
    // template).
    
    // Create instances of the class "UserParticleBoundaryFunctions" for different Ptypes

    // Q: Could we use THE_PARTICLE_TYPE here, instead of "cartesian_xy"?
    // A: Yes, if we were making a .so library for a specific Ptype only.
    
    // general version of .so
    makeUserParticleBoundaryFunctions<Ptype::PARTICLE_TYPE>(m, TOSTRING(PARTICLE_TYPE));
    
  } // ENDDEF: PYBIND11_MODULE()

} // namespace dnt
