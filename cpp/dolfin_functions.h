/*! \file dolfin_functions.h

  \brief This is the header file for dolfin_functions.cpp

  \namespace dnt
  \sa dolfin_functions.cpp

*/

#ifndef DOLFINFUNCTIONS_H
#define DOLFINFUNCTIONS_H

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
#include <pybind11/pybind11.h>

//#include <Eigen/Dense>

#include "Pstruct.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

namespace py = pybind11;

namespace dnt
{

  py::tuple find_facet(py::object mesh_M, double* x0, double* dx, size_t cell_index, bool returnStdArray = false);

  bool is_inside_vertices(dolfin::Mesh& mesh, const unsigned int* vertices, double* point);

  //template <Ptype PT, typename Ftype>
  template <Ptype PT>
  void interpolate_field_to_points(dolfin::Function* field,
                                   py::array_t<Pstruct<PT>, 0> points, // py::array_t
                                   py::ssize_t npoints,
                                   py::array_t<double> field_at_points);
  
  template <Ptype PT>
  void interpolate_field_to_points(dolfin::Function* field,
                                   Pstruct<PT>* points, // array of Pstructs
                                   py::ssize_t npoints,
                                   py::array_t<double> field_at_points);
  
} // namespace dnt

#endif
