// Copyright (C) 2018 L. D. Hughes
// pseg_cpp uses pybind11 to create Python-callable functions coded in C++ that
// perform computations on a single segment of the SegmentedArray_C data structure
// used to store particles in DnT.

// A SegmentedArray_C object is made up of Numpy record arrays called psegs.
// pybind11 allows us to directly access the data in a pseg.  The C++ functions speed
// up the computations on the pseg.

// Follow the example at https://github.com/pybind/pybind11/blob/534b756cb3244aca959359c07d44e4fed3498ba8/tests/test_numpy_dtypes.cpp

// Compiled on Xubuntu 18.04 with:

// c++ -O3 -Wall -shared -std=c++11 -fPIC -I/usr/include/eigen3 `python3 -m pybind11 --includes` pseg_cpp.cpp -o pseg_cpp`python3-config --extension-suffix`

// The first backquote above resolves to these include files: -I/usr/include/python3.6m -I/usr/local/include/python3.6 -I/home/tph/.local/include/python3.6m
// The second resolves to this extension: .cpython-36m-x86_64-linux-gnu.so
// 3.6m means the --with-pymalloc flag is on.

// Dolfin C++ source:
// GenericVector: /home/tph/workspace/dolfin/dolfin/la/

/*
Contents:
  Particle records in C++:
    DnT_prec1D struct and operator<<
    DnT_prec2D struct and operator<<
    DnT_prec3D struct and operator<<

  Print a pseg array:
    template <typename S>
    py::list print_recarray(py::array_t<S, 0> arr)

  Copy spatial coordinates from a particle record to a double array:
    template <typename S>
    void prec_to_double(S& p3D, double* point)

  Interpolate particle weights to DoFs:
    template <typename S>
    void interpolate_weights_to_dofs(py::array_t<S, 0> pseg, dolfin::Function& dF)

  Sum the particle weights in each cell:
    template <typename S>
    void add_weights_to_cells(py::array_t<S, 0> pseg, dolfin::Function& dF)

  PYBIND11 interfaces to Python for all the above 

*/

// #include <memory>

// #include <pybind11/pybind11.h>
// #include <pybind11/numpy.h>
// #include <pybind11/stl.h>
// #include <pybind11/eigen.h>
// #include <pybind11/operators.h>

#include "prec.h"

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

// Define the C++ function interpolate_weights_to_dofs(): it's templated on the type
// of the struct that corresponds to the Numpy pseg array. It loops over each
// particles in pseg and adds its contribution to the DoFs in the cell contraining
// the particle.
template <typename S>
void interpolate_weights_to_dofs(py::array_t<S, 0> pseg, dolfin::Function& dF) { // The 0 means a continguous array with C ordering

  const auto pseg_info = pseg.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
  const auto p = static_cast<S*>(pseg_info.ptr); // Pointer to a particle record in pseg
  
  std::shared_ptr< const dolfin::FunctionSpace > dFS = dF.function_space();
  
  // Variables for cell information
  const dolfin::Mesh& mesh = *dFS->mesh();
  std::vector<double> dof_coordinates;
  ufc::cell ufc_cell;

  // Variables for evaluating basis
  const std::size_t rank = dFS->element()->value_rank();
  std::size_t size_basis = 1;
  for (std::size_t i = 0; i < rank; ++i)
    size_basis *= dFS->element()->value_dimension(i);
  std::size_t dofs_per_cell = dFS->element()->space_dimension();
  std::vector<double> basis(size_basis);

// This vector holds the values contributed by a particle to the density array.
  std::vector<double> weights(dofs_per_cell);
  
  // Variables for adding local information to vector
  double basis_sum;
  
  std::shared_ptr< const dolfin::GenericDofMap > dDofmap = dFS->dofmap();

  // Print the DoF map. In this case it maps cell indices to DoF indices.
  // dFS->print_dofmap();

  // Add the particle weights to the Vector
  double point[] = {0.0, 0.0, 0.0};
  for (auto ip = 0; ip < pseg_info.size; ip++) {
    // std::cout << ip << " x=" << p[ip].x_ << std::endl;

    // Copy spatial coordinates of p[ip] to a double[] for passing to evaluate_basis()
    prec_to_double<S>(p[ip], point);
    
    auto cellIndex = p[ip].cell_index_;
    // std::cout << ip << " cellIndex=" << cellIndex << std::endl;

    // Get the indices of the DoFs in this cell.
    auto dofIndices = dDofmap->cell_dofs(cellIndex);
    
    // Get the cell object containing this particle
    dolfin::Cell cell(mesh, static_cast<std::size_t>(cellIndex));
    // Get the coordinates of the DoFs in this cell
    cell.get_coordinate_dofs(dof_coordinates);

    // Evaluate all basis functions at the point[]
    cell.get_cell_data(ufc_cell);
    for (std::size_t i = 0; i < dofs_per_cell; ++i)
    {
      dFS->element()->evaluate_basis(i, basis.data(),
                                     point,  // const double* x
                                     dof_coordinates.data(),
                                     ufc_cell.orientation);
      basis_sum = 0.0;
      for (const auto& v : basis)
        basis_sum += v;
      weights[i] = p[ip].weight_*basis_sum;
    }

    //dF.vector() is a std::shared_ptr<GenericVector>
    dF.vector()->add_local(weights.data(), dofs_per_cell, dofIndices.data());
    
  }
  
}

// Define the C++ function add_weights_to_cells(): it's templated on the type of the
// struct that corresponds to the Numpy pseg array. It loops over the particles in
// pseg and adds their weights to a cell-density array.
template <typename S>
void add_weights_to_cells(py::array_t<S, 0> pseg, dolfin::Function& dF) { // The 0 means a contiguous array with C ordering

  const auto pseg_info = pseg.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
  const auto p = static_cast<S*>(pseg_info.ptr); // Pointer to a particle record in pseg
  
  std::shared_ptr< const dolfin::FunctionSpace > dFS = dF.function_space();
  std::size_t dofs_per_cell = dFS->element()->space_dimension();

// This vector hold the values contributed by a particle to the density array. In the
// case of the cell density, a particle contributes to the one and only DoF in the
// cell, so we use the name 'weight' instead of 'weights':
  std::vector<double> weight(dofs_per_cell);
  
  std::shared_ptr< const dolfin::GenericDofMap > dDofmap = dFS->dofmap();

  // Print the DoF map. In this case it's cell indices to DoF indices.
  // dFS->print_dofmap();

  // Add the particle weights to the Vector
  for (auto i = 0; i < pseg_info.size; i++) {
    // std::cout << i << " x=" << p[i].x_ << std::endl;
    auto cellIndex = p[i].cell_index_;
    // std::cout << i << " cellIndex=" << cellIndex << std::endl;

    // In the case of a cell density, there's only one DoF in the cell, so we use the
    // name 'dofIndex' instead of 'dofIndices':
    auto dofIndex = dDofmap->cell_dofs(cellIndex);
    
    weight[0] = p[i].weight_;  // Only 1 DoF, so there's only 1 element in 'weight[]'.
    
    //dF.vector() is a std::shared_ptr<GenericVector>
    dF.vector()->add_local(weight.data(), dofs_per_cell, dofIndex.data());

// Q: would it be better to use the loop to fill arrays and then call this function
// only once, after the loop?

  }
  
}

// Below is the interface code that makes the above functions callable from Python

// The PYBIND11_MODULE() macro creates a function that will be called when an import
// statement is issued from within Python. The module name (pseg_cpp) is given as the
// first macro argument (it should not be in quotes). The second argument (m) defines
// a variable of type py::module which is the main interface for creating
// bindings. The method module::def() generates binding code that exposes the add()
// function to Python.

// // Create a variable 'm' of type py::module
PYBIND11_MODULE(pseg_cpp, m) {

// typeinfo may be registered before the dtype descriptor for scalar casts to work...

// Connect the Python symbol add_weights_to_cells() to the C++ function declared as
//       void add_weights_to_cells(py::array_t<S, 0> pseg, dolfin::Function& dF) {}
  m.def("add_weights_to_cells1D", &add_weights_to_cells<DnT_prec1D>);
  m.def("add_weights_to_cells2D", &add_weights_to_cells<DnT_prec2D>);
  m.def("add_weights_to_cells3D", &add_weights_to_cells<DnT_prec3D>);

// Connect the Python symbol interpolate_weights_to_dofs to the C++ function declared as
//       void interpolate_weights_to_dofs(py::array_t<S, 0> pseg, dolfin::Function& dF) {}
  m.def("interpolate_weights_to_dofs1D", &interpolate_weights_to_dofs<DnT_prec1D>);
  m.def("interpolate_weights_to_dofs2D", &interpolate_weights_to_dofs<DnT_prec2D>);
  m.def("interpolate_weights_to_dofs3D", &interpolate_weights_to_dofs<DnT_prec3D>);
  
}


