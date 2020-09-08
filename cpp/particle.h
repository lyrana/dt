/*! \file particle.h

  \brief This file has the source code for a C++ implementation of particle advance algorithms.

  \namespace dnt

  template <size_t N_CELL_FACETS>
    void advance_charged_species_in_E_field_cartesian_x()
    void advance_charged_species_in_E_field_cartesian_xy()

  template <size_t N_CELL_FACETS>    
    void advance_neutral_species_cartesian_xyz()

  template <Ptype PT>
    void advance_charged_species_in_uniform_fields_cartesian_xyz()

  template <Ptype PT>
    void makeParticleMeshBoundaryConditions(py::module &m, std::string const &PT_str)

  \sa particle_solib.cpp, MeshEntityArrays.h, SegmentedArrayPair.h, ParticleMeshBoundaryConditions.h, UserParticleBoundaryFunctions.h, user_particle_boundary_functions_solib.cpp
*/

#ifndef PARTICLE_H
#define PARTICLE_H

#include <dolfin/common/Array.h>
#include <dolfin/common/Hierarchical.h>
#include <dolfin/function/assign.h>
#include <dolfin/function/Constant.h>
#include <dolfin/function/Expression.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionAssigner.h>
#include <dolfin/function/FunctionAXPY.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/MultiMeshFunction.h>
#include <dolfin/function/MultiMeshFunctionSpace.h>
#include <dolfin/function/MultiMeshSubSpace.h>
#include <dolfin/function/LagrangeInterpolator.h>
#include <dolfin/function/SpecialFunctions.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/fem/MultiMeshDofMap.h>
#include <dolfin/geometry/Point.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/mesh/Mesh.h>

#include <iostream>
#include <tuple>
#include "SegmentedArrayPair.h"
#include "MeshEntityArrays.h"
#include "ParticleMeshBoundaryConditions.h"
#include "UserParticleBoundaryFunctions.h"

#include "dolfin_functions.h"

namespace py = pybind11;

// Static scratch space
static double particlePosition[3];
// A scratch array to hold: x,y,z, x0,y0,z0 (or a subset)
static double pCoord2[6];
// A scratch array to hold: dx,dy,dz (or a subset)
static double dx[3];

static dolfin::Cell* pcellPtr; // File scoped needed? Probably not.

// From https://en.cppreference.com/w/cpp/language/function_template:
// A function template by itself is not a type, or a function, or any other entity. No
// code is generated from a source file that contains only template definitions. In order
// for any code to appear, a template must be instantiated: the template arguments must be
// determined so that the compiler can generate an actual function (or class, from a class
// template).

namespace dnt
{

    //! Create objects needed to advance a particle species.
    /*!
      The Numpy arrays passed in have been allocated in Particle_Module.py.
      
      \return void

    */
  void initialize_particle_integration(py::array_t<double> negE_in, py::array_t<double> Eext_in, py::array_t<double> zeroE_in)
    {
      std::cout << "Hello from initialize_particle_integration@1" << std::endl;
      // Initialize the static interpolation arrays
    }
  // ENDDEF: void initialize_particle_integration

  // Charged particles

  //! Advance a cartesian_x charged-particle species by one time-increment on a mesh.
#include "advance_charged_species_in_E_field_cartesian_x.cpp"  
  //! Advance a spherical_r charged-particle species by one time-increment on a mesh.
#include "advance_charged_species_in_E_field_spherical_r.cpp"  

  //! Advance a cartesian_xy charged-particle species by one time-increment on a mesh.
#include "advance_charged_species_in_E_field_cartesian_xy.cpp"

  //! Advance a cartesian_xyz charged-particle species one timestep in uniform fields.
#include "advance_charged_species_in_uniform_fields_cartesian_xyz.cpp"

  // Neutral particles

  //! Advance a cartesian_xy neutral-particle species by one time-increment on a mesh.
#include "advance_neutral_species_cartesian_xy.cpp"
  
  //! Advance a cartesian_xyz neutral-particle species by one time-increment on a mesh.
#include "advance_neutral_species_cartesian_xyz.cpp"
  
} // namespace dnt


// Do these cause the compiler to make the specialized classes?
// The .def bindings in particle_solib.cpp seem to be sufficient to get compilation.
// That's probably because the definitions are included here in particle.h.

//template void dnt::advance_neutral_species<dnt::Ptype::cartesian_xyz, 2>(py::object,
//                                                                         py::str,
//                                                                         py::object);
//template void dnt::advance_neutral_species<dnt::Ptype::cartesian_xy, 3>(py::object,
//                                                                        py::str,
//                                                                        py::object);

// Write a helper function to make specialized ParticleMeshBoundaryConditions
// classes. This is in particle.h, reflecting the fact that the Python version is in
// Particle_Module.py. This also avoids putting a long piece of code into
// particle_solib.cpp, and avoids having to make a separate module for it.  Also,
// particle.h already has the pybind11 headers.

namespace dnt {

  // The anonymous namespace limits the scope of the functions in it to this file.
  namespace {

    //! Make ParticleMeshBoundaryConditions classes to contain UserParticleBoundaryFunctions for different Ptypes.
    /*!

      makeParticleMeshBoundaryConditions() is a 'helper' function. It creates the Python
      bindings for a specialized ParticleMeshBoundaryConditions class and its member
      functions based on the PT template parameter. Class instances can then be created
      and used in Python.

      \param PT is a template parameter for the function, specifying the class Ptype.
      \param m is a py::module object created by PYBIND11_MODULE.
      \param PT_str is a string used to create a unique Python class-name based on the Ptype.

      \return void

      \sa UserParticleBoundaryFunctions.h, UserParticleBoundaryFunctions.cpp

     */
    template <Ptype PT>
    void makeParticleMeshBoundaryConditions(py::module &m, std::string const &PT_str)
    {
      using PMBC = ParticleMeshBoundaryConditions<PT>;
      // Make a class name with the particle structure type (Ptype) string appended to the
      // string "ParticleMeshBoundaryConditions_"
      // std::string pyclass_name = std::string("ParticleMeshBoundaryConditions_") + PT_str;
      // Don't need the PT_str in the name if the MODULE_NAME has it.
      std::string pyclass_name = std::string("ParticleMeshBoundaryConditions");      
      
      // Create the Python binding for this class
      py::class_<PMBC>(m, pyclass_name.c_str())
        
        // The C++ ctor args types are template parameters in py::init<>()
        // Note that if a py::arg() specification is needed (e.g., to specify a
        // default value), then every argument has to have a py::arg().
        .def(py::init<std::vector<std::string>&, py::object&, UserParticleBoundaryFunctions<PT>&, bool>(), py::arg("species_names"), py::arg("pmesh_M"), py::arg("userParticleBoundaryFunctions"), py::arg("print_flag") = false);      
      // cf. ~/workspace/dolfin/python/src/geometry.cpp for __getitem__, if needed.
        
    } // void makeParticleMeshBoundaryConditions
    
  } // namespace is anonymous

}
#endif
