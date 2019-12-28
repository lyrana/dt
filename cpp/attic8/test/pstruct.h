// Copyright (C) 2018 L. D. Hughes

// A Numpy structured array corresponds to a struct in C/C++. The struct
// member-variable names have to be the same as the field names in the Python structured
// array unless PYBIND11_NUMPY_DTYPE_EX is used.  The following struct corresponds to
// a particle structured in Python.

/*
Contents:
  Particle structs in C++:
    DnT_pstruct1D struct and operator<<
    DnT_pstruct2D struct and operator<<
    DnT_pstruct3D struct and operator<<

  Put a pseg array in a Python list:
    template <typename PS>
    py::list print_pstructarray(py::array_t<PS, 0> arr)

  Copy spatial coordinates from a particle struct to a double array:
    template <typename PS>
    void pstruct_to_double(PS& ps, double* point)

  The PYBIND11 interfaces to Python for all the above.

*/

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>

namespace py = pybind11;

// For (x,) particle coordinates
struct DnT_pstruct1D {
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

// Declare the << operator for the DnT_pstruct1D data type. This is used in print_pstructarray()
std::ostream& operator<<(std::ostream& os, const DnT_pstruct1D& p);

// For (x, y) particle coordinates
struct DnT_pstruct2D {
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

// Declare the << operator for the DnT_pstruct2D data type. This is used in print_pstructarray()
std::ostream& operator<<(std::ostream& os, const DnT_pstruct2D& p);

// For (x, y, z) particle coordinates
struct DnT_pstruct3D {
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

// Declare the << operator for the DnT_pstruct3D data type. This is used in print_pstructarray()
std::ostream& operator<<(std::ostream& os, const DnT_pstruct3D& p);

// Declare the C++ function print_pstructarray(): it's templated on the type of the struct
// that corresponds to the Numpy structured array. It takes the array of structs, writes
// all the values into a string, and returns the string. This function is used later to
// define a Python-callable function to write the values in a structured array.
template <typename PS>
py::list print_pstructarray(py::array_t<PS, 0> arr);

// This function copies the spatial coordinates from a particle struct to a double array.
// It uses template specialization to handle particle structs with different dimensions.
template <typename PS>
void pstruct_to_double(PS& ps, double* point);
