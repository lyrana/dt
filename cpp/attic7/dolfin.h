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

#include <Eigen/Dense>

#include "Pstruct.h"
#include "Fstruct.h"

// Put these in namespace dnt:

namespace dnt
{
  
void divide_by_cell_volumes(dolfin::Function& dF, const std::map<int, double> &cell_volume_dict);
template <typename PS>
void interpolate_field_to_points(dolfin::Function& field,
                                 py::array_t<PS, 0> points,
                                 py::array_t<double> field_at_points);

/*
template <typename S>
void interpolate_field_to_points_v1(dolfin::Function& field,
                                    py::array_t<S, 0> points,
                                    py::array_t<double> field_at_points);
*/

//template <typename PS>
//bool cell_contains_point(dolfin::Mesh& mesh, py::array_t<int> vertices, PS point);

// Distinguish functions by their arguments
// These are NOT used!
//bool cell_contains_point(dolfin::Mesh& mesh, py::array_t<int> vertices, DnT_pstruct1D point);
//bool cell_contains_point(dolfin::Mesh& mesh, py::array_t<int> vertices, DnT_pstruct2D point);
 bool cell_contains_point(dolfin::Mesh& mesh, py::array_t<int> vertices, Pstruct<Ptype::cartesian_xyz> point);

// These ARE used!
bool cell_contains_point(dolfin::Mesh& mesh, py::array_t<int> vertices, double& x);
bool cell_contains_point(dolfin::Mesh& mesh, py::array_t<int> vertices, double& x, double& y);
bool cell_contains_point(dolfin::Mesh& mesh, py::array_t<int> vertices, double& x, double& y, double& z);

} // namespace dnt
