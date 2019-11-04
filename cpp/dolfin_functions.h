/*! \file dolfin_functions.h

  \brief This is the header file for dolfin_functions.cpp

  \namespace dnt
  \sa dolfin_functions.cpp

*/

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

//#include "Pstruct.h"
//#include "Fstruct.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
// #include <pybind11/eigen.h>
#include <pybind11/operators.h>
//#include <pybind11/cast.h>

namespace py = pybind11;


// Put these in namespace dnt:

namespace dnt
{
  
  bool cell_contains_point(dolfin::Mesh& mesh, const unsigned int* vertices, double* point);
  py::tuple find_facet(py::object mesh_M, double* x0, double* dx, size_t cell_index);
  
} // namespace dnt
