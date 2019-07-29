// Copyright (C) 2018 L. D. Hughes

// A Numpy record array corresponds to a struct in C/C++. The struct
// member-variable names have to be the same as the field names in the Python record
// array unless PYBIND11_NUMPY_DTYPE_EX is used.  The following struct corresponds to
// a particle record in Python.

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

  The PYBIND11 interfaces to Python for all the above.

*/

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>

/* #include <dolfin/common/Array.h> */
/* #include <dolfin/common/Hierarchical.h> */
/* #include <dolfin/function/assign.h> */
/* #include <dolfin/function/Constant.h> */
/* #include <dolfin/function/Expression.h> */
/* #include <dolfin/function/Function.h> */
/* #include <dolfin/function/FunctionAssigner.h> */
/* #include <dolfin/function/FunctionAXPY.h> */
/* #include <dolfin/function/FunctionSpace.h> */
/* #include <dolfin/function/MultiMeshFunction.h> */
/* #include <dolfin/function/MultiMeshFunctionSpace.h> */
/* #include <dolfin/function/MultiMeshSubSpace.h> */
/* #include <dolfin/function/LagrangeInterpolator.h> */
/* #include <dolfin/function/SpecialFunctions.h> */
/* #include <dolfin/fem/FiniteElement.h> */
/* #include <dolfin/fem/GenericDofMap.h> */
/* #include <dolfin/fem/MultiMeshDofMap.h> */
/* #include <dolfin/geometry/Point.h> */
/* #include <dolfin/la/GenericVector.h> */
/* #include <dolfin/la/GenericMatrix.h> */
/* #include <dolfin/mesh/Mesh.h> */

namespace py = pybind11;

// For (x,) particle coordinates
struct DnT_prec1D {
// Using PYBIND11_NUMPY_DTYPE_EX, the C++ variable names can be different from the Python names
  double x_;                  // f8    8 bytes
  double x0_;                 // f8   16 
  double ux_;                 // f8   24
  double weight_;             // f8   32
  int bitflags_;              // i4   36
  int cell_index_;            // i4   40
  int unique_ID_;             // i4   44
  int crossings_;             // i4   48
};

static const size_t PREC1D = 48;

// Define the << operator for the DnT_prec1D data type. This is used in print_recarray()
// below to print each element in the array.
std::ostream& operator<<(std::ostream& os, const DnT_prec1D& p) {
  return os <<   "x=" << p.x_
            << ", x0=" << p.x0_
            << ", ux=" << p.ux_
            << ", w=" << p.weight_
            << ", flgs=" << p.bitflags_
            << ", indx=" << p.cell_index_
            << ", id=" << p.unique_ID_
            << ", crossings=" << p.crossings_;
};

// For (x, y) particle coordinates
struct DnT_prec2D {
// Using PYBIND11_NUMPY_DTYPE_EX, the C++ variable names can be different from the Python names
  double x_;                  // f8    8 bytes
  double y_;                  // f8   16
  double x0_;                 // f8   24
  double y0_;                 // f8   32
  double ux_;                 // f8   40
  double uy_;                 // f8   48
  double weight_;             // f8   56
  int bitflags_;              // i4   60
  int cell_index_;            // i4   64
  int unique_ID_;             // i4   68
  int crossings_;             // i4   72
};
static const size_t PREC2D = 72;

// Define the << operator for the DnT_prec2D data type. This is used in print_recarray()
// below to print each element in the array.
std::ostream& operator<<(std::ostream& os, const DnT_prec2D& p) {
  return os <<   "x=" << p.x_
            << ", y=" << p.y_
            << ", x0=" << p.x0_
            << ", y0=" << p.y0_
            << ", ux=" << p.ux_
            << ", uy=" << p.uy_
            << ", w=" << p.weight_
            << ", flgs=" << p.bitflags_
            << ", indx=" << p.cell_index_
            << ", id=" << p.unique_ID_
            << ", crossings=" << p.crossings_;
};

// For (x, y, z) particle coordinates
struct DnT_prec3D {
// Using PYBIND11_NUMPY_DTYPE_EX, the C++ variable names can be different from the Python names
  double x_;                  // f8    8 bytes
  double y_;                  // f8   16
  double z_;                  // f8   24
  double x0_;                 // f8   32
  double y0_;                 // f8   40
  double z0_;                 // f8   48
  double ux_;                 // f8   56
  double uy_;                 // f8   64
  double uz_;                 // f8   72
  double weight_;             // f8   80
  int bitflags_;              // i4   84
  int cell_index_;            // i4   88
  int unique_ID_;             // i4   92
  int crossings_;             // i4   96
};
static const size_t PREC3D = 96;

// Define the << operator for the DnT_prec3D data type. This is used in print_recarray()
// below to print each element in the array.
std::ostream& operator<<(std::ostream& os, const DnT_prec3D& p) {
  return os <<   "x=" << p.x_
            << ", y=" << p.y_
            << ", z=" << p.z_
            << ", x0=" << p.x0_
            << ", y0=" << p.y0_
            << ", z0=" << p.z0_
            << ", ux=" << p.ux_
            << ", uy=" << p.uy_
            << ", uz=" << p.uz_
            << ", w=" << p.weight_
            << ", flgs=" << p.bitflags_
            << ", indx=" << p.cell_index_
            << ", id=" << p.unique_ID_
            << ", crossings=" << p.crossings_;
};

// Define the C++ function print_recarray(): it's templated on the type of the struct
// that corresponds to the Numpy record array. It takes the array of structs, writes
// all the values into a string, and returns the string. This function is used later to
// define a Python-callable function to write the values in a record array.
template <typename S>
py::list print_recarray(py::array_t<S, 0> arr) {
  const auto arr_info = arr.request();  // request() returns metadata about the array (ptr, ndim, size, shape)
  const auto p = static_cast<S*>(arr_info.ptr); // Pointer to a particle record in pseg

// See ~/local/include/pybind11/buffer_info.h
// arr_info attributes are:
// itemsize, size, format, ndim, shape, strides 
  std::cout << "  itemsize=" << arr_info.itemsize
            << ", size=" << arr_info.size
            << ", format=" << arr_info.format
            << ", ndim=" << arr_info.ndim
            << std::endl;

  std::cout << "The size the pseg record is " << sizeof(S) << " bytes" << std::endl;
  
  auto l = py::list();
  for (ssize_t i = 0; i < arr_info.size; i++) {
    std::stringstream ss;
    ss << p[i];  // Use the above definition of "<<" to write the i'th record as a string.
    l.append(py::str(ss.str()));
  }
  return l;
}

// This function copies the spatial coordinates from a particle record to a double array.
// It uses template specialization to handle particle records with different dimensions.
template <typename S>
void prec_to_double(S& p3D, double* point)
{
  // There is no generic version
}
template <>
void prec_to_double<DnT_prec1D>(DnT_prec1D& p1D, double* point)
{
  point[0] = p1D.x_;
}
template <>
void prec_to_double<DnT_prec2D>(DnT_prec2D& p2D, double* point)
{
  point[0] = p2D.x_;
  point[1] = p2D.y_;
}
template <>
void prec_to_double<DnT_prec3D>(DnT_prec3D& p3D, double* point)
{
  point[0] = p3D.x_;
  point[1] = p3D.y_;
  point[2] = p3D.z_;
}

// Below is the interface code that makes the above functions callable from Python

// The PYBIND11_MODULE() macro creates a function that will be called when an import
// statement is issued from within Python. The module name (pseg_cpp) is given as the
// first macro argument (it should not be in quotes). The second argument (m) defines
// a variable of type py::module which is the main interface for creating
// bindings. The method module::def() generates binding code that exposes the C++
// functions to Python.

// // Create a variable 'm' of type py::module
PYBIND11_MODULE(prec_cpp, m) {

// typeinfo may be registered before the dtype descriptor for scalar casts to work...

//  py::class_() creates Python bindings for a C++ class or struct-style data structure
  py::class_<DnT_prec1D>(m, "DnT_prec1D");
  py::class_<DnT_prec2D>(m, "DnT_prec2D");
  py::class_<DnT_prec3D>(m, "DnT_prec3D");

// Register DnT_prec1D as a Numpy dtype descriptor. DnT_prec1D is of type py::dtype (a class).
//  PYBIND11_NUMPY_DTYPE(DnT_prec1D, x, x0, ux, weight, bitflags, cell_index, unique_ID, crossings);
// The "_EX" variation allows the Python names of the variables in the record to be different from the variable names in the C++ struct.
  PYBIND11_NUMPY_DTYPE_EX(DnT_prec1D, x_, "x", x0_, "x0", ux_, "ux", weight_, "weight", bitflags_, "bitflags", cell_index_, "cell_index", unique_ID_, "unique_ID", crossings_, "crossings");

// Register DnT_prec2D and DnT_prec3D
  PYBIND11_NUMPY_DTYPE_EX(DnT_prec2D, x_, "x", y_, "y", x0_, "x0", y0_, "y0", ux_, "ux", uy_, "uy", weight_, "weight", bitflags_, "bitflags", cell_index_, "cell_index", unique_ID_, "unique_ID", crossings_, "crossings");

  PYBIND11_NUMPY_DTYPE_EX(DnT_prec3D, x_, "x", y_, "y", z_, "z", x0_, "x0", y0_, "y0", z0_, "z0", ux_, "ux", uy_, "uy", uz_, "uz", weight_, "weight", bitflags_, "bitflags", cell_index_, "cell_index", unique_ID_, "unique_ID", crossings_, "crossings");
  
/* now DnTprec1D, etc., can be used as template arguments to py::array_t */

// Connect the Python symbol print_pseg() to the C++ function declared as:
//                   template <typename S>
//                   py::list print_recarray(py::array_t<S, 0> arr) {}
  m.def("print_pseg1D", &print_recarray<DnT_prec1D>);
  m.def("print_pseg2D", &print_recarray<DnT_prec2D>);
  m.def("print_pseg3D", &print_recarray<DnT_prec3D>);

// Example: Connect the Python symbol f_simple() to the C++ function in braces {...}
// This prints only one record of the array
  m.def("f_simple", [](DnT_prec1D p) { return p.x_ * 10; });

}


