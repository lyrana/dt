// Copyright (C) 2019 L. D. Hughes

// Pass Python class objects to C++ and extract data from them.

/*
Contents:

  void function_with_several_args(bool tf, py::list& list3, py::array_t<double>& darray3)
  void function_with_DTcontrol_C_arg(py::object ctrl)
  void function_with_Particle_C_arg(py::object particle_P)
  void function_with_pseg_arg(py::ssize_t np_seg, py::array_t<Pstruct<PT>,0> pseg)

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
    void pstruct_to_point(PS& ps, double* point)

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

#include "SegmentedArrayPair.h"

namespace py = pybind11;

namespace dnt
{
  //! Test passing of various Python types and classes to C++

  void function_with_several_args(bool tf, py::list& int_list3, py::list& string_list3, py::array_t<double>& darray3)
  {
    std::cout << "Hello from C++: function_with_several_args" << std::endl;
    std::cout << "tf is " << tf << std::endl;
    
    for (size_t i = 0; i < int_list3.size(); i++)
      {
        std::cout << "int_list3[" << i << "] = " << int_list3[i].cast<int>() << std::endl;        
      }

    for (size_t i = 0; i < string_list3.size(); i++)
      {
        std::cout << "string_list3[" << i << "] = " << string_list3[i].cast<std::string>() << std::endl;
        // Note: std::string(py::str(string_list3[i])) also works.
      }

    auto darray3Proxy = darray3.unchecked<1>();
    // If we use size_t i, get warning: comparison between signed and unsigned integer expressions
    // py::ssize_t is signed.
    for (py::ssize_t i = 0; i < darray3Proxy.shape(0); i++)
      {
        
        std::cout << "darray3[" << i << "] = " << darray3Proxy(i) << std::endl;
      }

  }
  
  void function_with_DTcontrol_C_arg(py::object ctrl)
  {
    std::cout << "Hello from C++: function_with_DTcontrol_C_arg" << std::endl;
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
  // ENDDEF: void function_with_DTcontrol_C_arg(py::object ctrl)
  
  void function_with_Particle_C_arg(py::object particle_P)
  {
    std::cout << "Hello from C++: function_with_Particle_C_arg" << std::endl;

    // Access particle_P class variables:
    auto NO_FACET = particle_P.attr("pmesh_M").attr("NO_FACET").cast<int>();
    auto msg = "The value of NO_FACET is " + std::to_string(NO_FACET);
    std::cout << msg << std::endl;

    // Access the dolfin::Mesh object from the particle_P object:
    auto pmesh_M = particle_P.attr("pmesh_M");
    auto mesh = pmesh_M.attr("mesh").cast<dolfin::Mesh>();

    auto pDim = particle_P.attr("particle_dimension").cast<int>();
    msg = "pDim = " + std::to_string(pDim);
    std::cout << msg << std::endl;
    
    // Get the list of species names from particle_P as a py::list
    // Is the cast needed? NO, but then we can't use 'auto' as the type: it must be explicit.
    //    auto neutralSpecies = particle_P.attr("neutral_species").cast<py::list>();
    py::list neutralSpecies = particle_P.attr("neutral_species");
    // Access the name of a species from the neutralSpecies py::list. Note the cast
    // is needed to get an std::string
    msg = "neutral species 0 = " + neutralSpecies[0].cast<std::string>();
    std::cout << msg << std::endl;

    // Get the list of species names from particle_P as an std::vector of strings
    auto neutralSpecies2 = particle_P.attr("neutral_species").cast<std::vector<std::string>>();
    // Access the name of a species from the neutralSpecies2 std::vector:
    msg = "neutral species 0 #2 = " + neutralSpecies2[0];
    std::cout << msg << std::endl;
     
    // Addess the SAPs in particle_P by casting to a py::dict
    
    auto sapDict = particle_P.attr("sap_dict").cast<py::dict>();

    // Make a loop over all the species in sapDict
    msg = "Loop on species:";
    std::cout << msg << std::endl;
    for (auto species : sapDict)
      {
        auto msg = "\tThis is species " + std::string(py::str(species.first));
        std::cout << msg << std::endl;
      }

    // Access a particular SAP in sapDict
    auto sapSp0 = sapDict[neutralSpecies[0]].cast<SegmentedArrayPair<Ptype::cartesian_xyz> *>();
    // Call a SAP function
    msg = "Length of segments in sapSp0 = " + std::to_string(sapSp0->get_segment_length());
    std::cout << msg << std::endl;
    
    // Access the SAPs in particle_P by casting to a std::map of species names
    // (std::strings) to references to a SAP of the right type.
    
    // Note the pointer star: .cast<MyClass *>, not .cast<Myclass>. This is because
    // sapDict is a Python dictionary and the dictionary values are *references* to SAPs.
    auto sapMap = particle_P.attr("sap_dict").cast<std::map<std::string, SegmentedArrayPair<Ptype::cartesian_xyz> *>>();

    auto sapH = sapMap["neutral_H"];
    msg = "Length of segments in sapH = " + std::to_string(sapH->get_segment_length());
    std::cout << msg << std::endl;
    
  }
  // ENDDEF: void function_with_Particle_C_arg(py::object particle_P)

  
  template <Ptype PT>
    void function_with_pseg_arg(py::ssize_t np_seg, py::array_t<Pstruct<PT>,0> pseg)
  {
    std::cout << "Hello from C++: function_with_pseg_arg" << std::endl;
    std::cout << "There are " << np_seg << " elements in the passed pseg array" << std::endl;

    auto psegProxy = pseg.unchecked();
    for (auto ip = 0; ip < np_seg; ip++) {
      std::cout << "Particle " << ip << std::endl;
      std::cout << "\tx=" << psegProxy(ip).x_ << std::endl;
      std::cout << "\ty=" << psegProxy(ip).y_ << std::endl;
      std::cout << "\tz=" << psegProxy(ip).z_ << std::endl;
      std::cout << "\tx0=" << psegProxy(ip).x0_ << std::endl;
      std::cout << "\ty0=" << psegProxy(ip).y0_ << std::endl;
      std::cout << "\tz0=" << psegProxy(ip).z0_ << std::endl;
      std::cout << "\tux=" << psegProxy(ip).ux_ << std::endl;
      std::cout << "\tuy=" << psegProxy(ip).uy_ << std::endl;
      std::cout << "\tuz=" << psegProxy(ip).uz_ << std::endl;
      std::cout << "\tweight=" << psegProxy(ip).weight_ << std::endl;
      std::cout << "\tcell_index=" << psegProxy(ip).cell_index_ << std::endl;
      std::cout << "\tcrossings=" << psegProxy(ip).crossings_ << std::endl;
    }

  }
  
} // namespace dnt

#endif
