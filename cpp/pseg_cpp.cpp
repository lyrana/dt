// pseg_cpp uses pybind11 to create Python-callable functions coded in C++ that perform computations
// on a single segment of the SegmentedArray_C data structure used to store particles in DnT.

// A SegmentedArray_C object is made up of Numpy record arrays called psegs.
// pybind11 allows us to directly access the data in a pseg.  The C++ functions speed
// up the computations on the pseg.

// Follow the example at https://github.com/pybind/pybind11/blob/534b756cb3244aca959359c07d44e4fed3498ba8/tests/test_numpy_dtypes.cpp

// Compiled on Xubuntu 18.04 with:

// c++ -O3 -Wall -shared -std=c++11 -fPIC -I/usr/include/eigen3 `python3 -m pybind11 --includes` pseg_cpp.cpp -o pseg_cpp`python3-config --extension-suffix`

// The first backtick above resolves to these include files: -I/usr/include/python3.6m -I/usr/local/include/python3.6 -I/home/tph/.local/include/python3.6m
// The second resolves to this extension: .cpython-36m-x86_64-linux-gnu.so
// 3.6m means the --with-pymalloc flag is on.

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

namespace py = pybind11;

// A Numpy record array corresponds to a struct in C/C++. The struct
// member-variable names have to be the same as the field names in the Python record
// array unless PYBIND11_NUMPY_DTYPE_EX is used.  The following struct corresponds to
// a particle record in Python.
struct DnT_prec {
  
  // double x;
  // double x0;
  // double ux;
  // double weight;
  // int bitflags;
  // int cell_index;
  // int unique_ID;
  // int crossings;

// Using PYBIND11_NUMPY_DTYPE_EX, the member variable names can be different from the Python names  
  double x_;
  double x0_;
  double ux_;
  double weight_;
  int bitflags_;
  int cell_index_;
  int unique_ID_;
  int crossings_;
};

// Define the << operator for the DnTprec data type. This is used in print_recarray()
// below to print each element in the array.
std::ostream& operator<<(std::ostream& os, const DnT_prec& p) {
//  return os << "p:" << p.x << "," << p.x0 << "," << p.weight << "," << p.bitflags << "," << p.cell_index << "," << p.unique_ID << "," << p.crossings;
  return os << "p:" << p.x_ << "," << p.x0_ << "," << p.weight_ << "," << p.bitflags_ << "," << p.cell_index_ << "," << p.unique_ID_ << "," << p.crossings_;
}

// Define the C++ function print_recarray(): it's templated on the type of the struct
// that corresponds to the Numpy record array. It takes the array of structs, writes
// all the values into a string, and returns the string. This function is used later to
// define a Python-callable function to write the values in a record array.
template <typename S>
py::list print_recarray(py::array_t<S, 0> arr) {
  const auto req = arr.request();
  const auto ptr = static_cast<S*>(req.ptr);
  auto l = py::list();
  for (ssize_t i = 0; i < req.size; i++) {
    std::stringstream ss;
    ss << ptr[i];
    l.append(py::str(ss.str()));
  }
  return l;
}

// Define the C++ function add_weights_to_cells(): it's templated on the type of the
// struct that corresponds to the Numpy record array. It ...
template <typename S>
void add_weights_to_cells(py::array_t<S, 0> pseg, dolfin::Function& dF) { // The 0 means a continguous array with C ordering
  const auto pseg_info = pseg.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
  const auto ptr = static_cast<S*>(pseg_info.ptr);

//  auto l = py::list();

// Loop on the particles
//  for (size_t i = 0; i < pseg_info.size; i++) {
  for (auto i = 0; i < pseg_info.size; i++) {
    std::cout << i << "x" << ptr[i].x_;
  }
}


// Below is the code that makes the above functions callable from Python

// The PYBIND11_MODULE() macro creates a function that will be called when an import
// statement is issued from within Python. The module name (pseg_cpp) is given as the
// first macro argument (it should not be in quotes). The second argument (m) defines
// a variable of type py::module which is the main interface for creating
// bindings. The method module::def() generates binding code that exposes the add()
// function to Python.

// // Create a variable 'm' of type py::module
PYBIND11_MODULE(pseg_cpp, m) {

// typeinfo may be registered before the dtype descriptor for scalar casts to work...

//  py::class_() creates bindings for a C++ class or struct-style data structure
  py::class_<DnT_prec>(m, "DnT_prec");

// Register DnT_prec as a Numpy dtype descriptor. DnT_prec is of type py::dtype (a class).
//  PYBIND11_NUMPY_DTYPE(DnT_prec, x, x0, ux, weight, bitflags, cell_index, unique_ID, crossings);
// The "_EX" variation allows the Python names of the variables in the record to be different from the variable names in the C++ struct.
  PYBIND11_NUMPY_DTYPE_EX(DnT_prec, x_, "x", x0_, "x0", ux_, "ux", weight_, "weight", bitflags_, "bitflags", cell_index_, "cell_index", unique_ID_, "unique_ID", crossings_, "crossings");

/* now both A and B can be used as template arguments to py::array_t */

// Make Python symbol print_pseg() call C++ function print_recarray()
  m.def("print_pseg", &print_recarray<DnT_prec>);

// Make Python function f_simple() call C++ function {...}.
// This prints only one record
//  m.def("f_simple", [](DnT_prec p) { return p.x * 10; });
  m.def("f_simple", [](DnT_prec p) { return p.x_ * 10; });

  m.def("add_weights_to_cells", &add_weights_to_cells<DnT_prec>);
  
}


