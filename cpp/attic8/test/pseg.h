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

#include "pstruct.h"

template <typename PS>
void add_weights_to_cells(py::array_t<PS, 0> pseg, dolfin::Function& dF);

template <typename PS>
void interpolate_weights_to_dofs(py::array_t<PS, 0> pseg, dolfin::Function& dF);

//template <>
//void add_weights_to_cells<DnT_pstruct1D>(py::array_t<DnT_pstruct1D, 0> pseg, dolfin::Function& dF);

//void interpolate_weights_to_dofs(py::array_t<DnT_pstruct1D, 0> pseg, dolfin::Function& dF);
