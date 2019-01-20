// Copyright (C) 2018 L. D. Hughes
// dolfin_cpp uses pybind11 to create Python-callable functions coded in C++ that
// perform computations on Dolfin objects used in DnT.

// Compiled on Xubuntu 18.04 with:

// c++ -O3 -Wall -shared -std=c++11 -fPIC -I/usr/include/eigen3 `python3 -m pybind11 --includes` dolfin_cpp.cpp -o dolfin_cpp`python3-config --extension-suffix`

// The first backquote above resolves to these include files: -I/usr/include/python3.6m -I/usr/local/include/python3.6 -I/home/tph/.local/include/python3.6m
// The second resolves to this extension: .cpython-36m-x86_64-linux-gnu.so
// 3.6m means the --with-pymalloc flag is on.

// Dolfin C++ source:
// GenericVector: /home/tph/workspace/dolfin/dolfin/la/

/*
Contents:

  Print a Python dictionary:
    void print_dict(py::dict dict)

  Divide the values of a cell function by the cell volumes:
    void divide_by_cell_volumes(dolfin::Function& dF, const std::map<int, double> &cell_volume_dict)

*/

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>

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

#include <typeinfo>

namespace py = pybind11;

//! Print the key:value pairs in a Python dictionary
void print_dict(py::dict dict) {
    for (auto item : dict)
    {
        std::cout << "key=" << std::string(py::str(item.first)) << ", "
                  << "value=" << std::string(py::str(item.second)) << std::endl;

// These give N8pybind116handleE as the type.
//        std::cout << "type of item.first: " << typeid(item.first).name() << std::endl;
//        std::cout << "type of item.second: " << typeid(item.second).name() << std::endl;
    }

}

//! Divide the values of a cell function by the cell volumes.
/*!
    \param dF is a dolfin function defined on the cells of a mesh. dF.vector() is a
    std::shared_ptr<GenericVector>.

    \param cell_volume_dict is a reference to a python dict containing the cell volumes.

    \return void
*/
void divide_by_cell_volumes(dolfin::Function& dF, const std::map<int, double> &cell_volume_dict) {
  
    std::shared_ptr< const dolfin::FunctionSpace > dFS = dF.function_space();
    std::size_t dofs_per_cell = dFS->element()->space_dimension();

    // In the case of the cell density, there's only one DoF in the cell, so we use
    // the name 'value' instead of 'values':

    std::vector<double> value(dofs_per_cell);
  
    std::shared_ptr< const dolfin::GenericDofMap > dDofmap = dFS->dofmap();
      
    auto mesh = dFS->mesh();
    // Loop on cells
    /// A CellIterator is a MeshEntityIterator of topological codimension 0.
    /// typedef MeshEntityIteratorBase<Cell> CellIterator;
    for (dolfin::CellIterator cell(*mesh); !cell.end(); ++cell)
    {
      // Get the dofIndex of the cell
      auto cellIndex = cell->index();
      // If there's only one DoF in the cell, so we use the
      // name 'dofIndex' instead of 'dofIndices':
      auto dofIndex = dDofmap->cell_dofs(cellIndex);
      dF.vector()->get_local(value.data(), dofs_per_cell, dofIndex.data());
      
      // Get the volume of the cell
      auto cellVol = cell_volume_dict.at(cellIndex);

      // This gives type "d" (double?)
      // std::cout << "type of cellVol: " << typeid(cellVol).name() << std::endl;
      
      // Divide the value by the volume
      value[0] /= cellVol;
      // Put back the value
      dF.vector()->set_local(value.data(), dofs_per_cell, dofIndex.data());        
    }
}

//! Determine if a point is inside a cell
/*!
    \param dF is a dolfin function defined on the cells of a mesh. dF.vector() is a
    std::shared_ptr<GenericVector>.

    \param cell_vertex_dict is a reference to a python dict containing cell vertices.

    \return void
*/

// void cell_contains_point(const std::map<int, dolfin::Cell&> &cell_dict) {}

/*
void cell_contains_point(dolfin::Mesh& mesh, const std::map<std::size_t, const dolfin::Cell& c> &cell_dict, std::size_t icell, double point)
{
  // Get the vertices of this cell
  auto cell = cell_volume_dict[icell];
  
  p0 = 0.0;
  mesh.coordinates();

}
*/

 //bool cell_contains_point_1d(dolfin::Mesh& mesh, py::list& vertices, double point)
 //bool cell_contains_point_1d(dolfin::Mesh& mesh, py::array_t<dolfin::uint32> vertices, double point)

 //bool cell_contains_point_1d(dolfin::Mesh& mesh, py::array_t<unsigned int> vertices, double point)

bool cell_contains_point_1d(dolfin::Mesh& mesh, unsigned int v0, unsigned int v1, double point)
{
//  std::size_t v0 = vertices[0];
//  std::size_t v1 = vertices[1];
  
// Use this for a Python list of ints called vertices[]:
  
//  auto v0 = vertices[0].cast<int>();
//  auto v1 = vertices[1].cast<int>();
  
//  auto v0 = vertices[0];
//  auto v1 = vertices[1];

//  cout << "v0: " << v0 << "v1: " << v1 << std::endl;
  
  double p0 = mesh.coordinates()[v0];
  double p1 = mesh.coordinates()[v1];
  
  if (p0 > p1)
    std::swap(p0, p1);
  return p0 <= point and point <= p1;
}

bool cell_contains_point_1d_bak(double p0, double p1, double point)
{
  if (p0 > p1)
    std::swap(p0, p1);
  return p0 <= point and point <= p1;
}


// Evaluate a dolfin Function at a set of points      
// This is the C++ version of a Python function in Dolfin_Module.py:
//     interpolate_field_to_points(self, points, field_at_points)
template <typename S>
void interpolate_field_to_points(dolfin::Function& field,
                                 py::array_t<S, 0> points,
                                 py::array_t<double> field_at_points)
{
  const auto points_info = points.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
  const auto p = static_cast<S*>(points_info.ptr); // Pointer to a record in points
  
  std::shared_ptr< const dolfin::FunctionSpace > fieldSpace = field.function_space();
  
  // Variables for cell information
  const dolfin::Mesh& mesh = *fieldSpace->mesh();
//  std::vector<double> dof_coordinates;
  ufc::cell ufc_cell;

  double point[] = {0.0, 0.0, 0.0};
  for (auto ip = 0; ip < points_info.size; ip++) {
    // std::cout << ip << " x=" << p[ip].x_ << std::endl;

    // Copy spatial coordinates of p[ip] to a double[] for passing to evaluate_basis()
    // ?? Can p[ip] be passed directly?
    // Or write a special form of eval!
    prec_to_double<S>(p[ip], point);
    
    auto cellIndex = p[ip].cell_index_;
    // std::cout << ip << " cellIndex=" << cellIndex << std::endl;

    // Get the indices of the DoFs in this cell.
//    auto dofIndices = dDofmap->cell_dofs(cellIndex);
    
    // Get the cell object containing this particle
    dolfin::Cell cell(mesh, static_cast<std::size_t>(cellIndex));
    cell.get_cell_data(ufc_cell);

    // Get the coordinates of the DoFs in this cell
    // cell.get_coordinate_dofs(dof_coordinates);

    // Evaluate the function at point
    // See dolfin/dolfin/function/Function.cpp
    field::eval(values, point, cell, ufc_cell);

/*    
    for (std::size_t i = 0; i < dofs_per_cell; ++i)
    {
      fieldSpace->element()->evaluate_basis(i, basis.data(),
                                     point,  // const double* x
                                     dof_coordinates.data(),
                                     ufc_cell.orientation);
      basis_sum = 0.0;
      for (const auto& v : basis)
        basis_sum += v;
      weights[i] = p[ip].weight_*basis_sum;
    }

    //field.vector() is a std::shared_ptr<GenericVector>
    field.vector()->add_local(weights.data(), dofs_per_cell, dofIndices.data());
*/
    
  }

}




// Below is the interface code that makes the above functions callable from Python

// The PYBIND11_MODULE() macro creates a function that will be called when an import
// statement is issued from within Python. The module name (dolfin_cpp) is given as the
// first macro argument (it should not be in quotes). The second argument (m) defines
// a variable of type py::module which is the main interface for creating
// bindings. The method module::def() generates binding code that exposes the C++
// functions to Python.

// Create a variable 'm' of type py::module
PYBIND11_MODULE(dolfin_cpp, m) {

// Connect the Python symbol print_dict() to the C++ function declared as
//         void print_dict(py::dict dict)  
  m.def("print_dict", &print_dict);
  
// Connect the Python symbol divide_by_cell_volumes() to the C++ function declared as
//         void divide_by_cell_volumes(dolfin::Function& dF, py::dict dict)
  m.def("divide_by_cell_volumes", &divide_by_cell_volumes);

  m.def("cell_contains_point_1d", &cell_contains_point_1d);
}


