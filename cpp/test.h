// Copyright (C) 2019 L. D. Hughes

// Tests for pybind11 in C++

/*
Contents:

  void function_with_DTcontrol_arg(py::object ctrl)
  void function_with_particle_P_arg(py::object particle_P)

//remove:
  Particle structs in C++:
    dnt::pstruct1D struct and operator<<
    dnt::pstruct2D struct and operator<<
    dnt::pstruct3D struct and operator<<

  Put a pseg array in a Python list:
    template <typename PS>
    py::list print_pstructarray(py::array_t<PS, 0> arr)

  Copy spatial coordinates from a particle struct to a double array:
    template <typename PS>
    void pstruct_to_double(PS& ps, double* point)

*/

#ifndef TEST_H
#define TEST_H

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


namespace py = pybind11;

namespace dnt
{
  //! Test passing of various DnT Python classes to C++
  
  void function_with_DTcontrol_arg(py::object ctrl)
  {
    std::cout << "Hello from function_with_DTcontrol_arg" << std::endl;
    auto title = ctrl.attr("title").cast<std::string>();
    auto author = ctrl.attr("author").cast<std::string>();
    auto dt = ctrl.attr("dt").cast<double>();
    auto dt_str = std::to_string(dt);
    auto n_timesteps = ctrl.attr("n_timesteps").cast<int>();
    auto n_timesteps_str = std::to_string(n_timesteps);
    
    auto msg = "The title is: " + title;
    std::cout << msg << std::endl;

    msg = "The author is: " + author;
    std::cout << msg << std::endl;

    msg = "The timestep is: " + dt_str;
    std::cout << msg << std::endl;

    msg = "The number of timesteps is: " + n_timesteps_str;
    std::cout << msg << std::endl;
  }

  void function_with_particle_P_arg(py::object particle_P)
  {
    std::cout << "Hello from function_with_particle_P_arg" << std::endl;

    //dolfin::Mesh& mesh
    
    auto NO_FACET = particle_P.attr("pmesh_M").attr("NO_FACET").cast<int>();

    auto pmesh_M = particle_P.attr("pmesh_M");
    auto pDim = particle_P.attr("particle_dimension").cast<int>();

    // Is the cast needed? NO, but then we can't use 'auto' as the type: it must be explicit.
    //    auto neutralSpecies = particle_P.attr("neutral_species").cast<py::list>();
    py::list neutralSpecies = particle_P.attr("neutral_species");

    //std::vector<int> &v
    auto neutralSpecies2 = particle_P.attr("neutral_species").cast<std::vector<std::string>>();
    
    /*    
void print_dict(py::dict dict) {
for (auto item : dict)
std::cout << "key=" << std::string(py::str(item.first)) << ", "
<< "value=" << std::string(py::str(item.second)) << std::endl;
}
    */
    auto pseg_arr = particle_P.attr("pseg_arr").cast<py::dict>();

    // Make a loop over all the species in pseg_arr
    for (auto species : pseg_arr)
      {
        auto msg = "This is species " + std::string(py::str(species.first));
        std::cout << msg << std::endl;
      }

    //    std::map<std::string, SegmentedArrayPair<Ptype::cartesian_x_y_z>> pseg_arr_map = particle_P.attr("pseg_arr").cast<std::map<std::string, SegmentedArrayPair<Ptype::cartesian_x_y_z>>>();
    // pseg.arr is cast to a std::map from species names (strings) to an SAP of the right type.
    // Note pointer: .cast<MyClass *>, not .cast<Myclass> !
    auto pseg_arr_map = particle_P.attr("pseg_arr").cast<std::map<std::string, SegmentedArrayPair<Ptype::cartesian_x_y_z> *>>();

    // compiler seems to turn this into a call to the SAP ctor.
    // Maybe that's what .cast<>() does: copies?
    // Fix: Use pointer-to-type above instead of type; see pybind11 docs p. 99.
    auto psa1 = pseg_arr_map["neutral_H"];
    auto msg = "Length of segments in psa1 = " + std::to_string(psa1->get_segment_length());
    std::cout << msg << std::endl;
    
    // Access the dolfin::Mesh object from the particle_P object:
    auto mesh = pmesh_M.attr("mesh").cast<dolfin::Mesh>();

    // Access particle_P class variables:
    msg = "The value of NO_FACET is " + std::to_string(NO_FACET);
    std::cout << msg << std::endl;
    
    msg = "pDim = " + std::to_string(pDim);
    std::cout << msg << std::endl;

    // Access the name of a species from the neutralSpecies py::list:
    msg = "neutral species 0 = " + neutralSpecies[0].cast<std::string>();
    std::cout << msg << std::endl;

    // Access the name of a species from the neutralSpecies2 std::vector:
    msg = "neutral species 0 #2 = " + neutralSpecies2[0];
    std::cout << msg << std::endl;

    // Make a loop over all the species:
    for (auto species : pseg_arr)
      {
        auto msg = "This is species " + std::string(py::str(species.first));
        std::cout << msg << std::endl;
      }

    // Access the SAP from py::ojbect pseg_arr using a cast
    auto psa = pseg_arr[neutralSpecies[0]].cast<SegmentedArrayPair<Ptype::cartesian_x_y_z> *>();
    // Call an SAP function
    msg = "Length of segments in psa = " + std::to_string(psa->get_segment_length());
    std::cout << msg << std::endl;
      
    /*    
    auto dt = ctrl.attr("dt").cast<double>();
    auto dt_str = std::to_string(dt);
    auto n_timesteps = ctrl.attr("n_timesteps").cast<int>();
    auto n_timesteps_str = std::to_string(n_timesteps);
    
    auto msg = "The author is " + author;
    std::cout << msg << std::endl;

    msg = "The timestep is " + dt_str;
    std::cout << msg << std::endl;

    msg = "The number of timesteps is " + n_timesteps_str;
    std::cout << msg << std::endl;
    */
  }

} // namespace dnt

#endif
