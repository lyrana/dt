// Copyright (C) 2018, 2019 L. D. Hughes

// A Numpy structured array corresponds to an array of struct in C/C++. The struct
// member-variable names have to be the same as the field names in the Python structured
// array unless PYBIND11_NUMPY_DTYPE_EX is used.  The following struct corresponds to
// a particle structure in Python.

/*
Contents:
  Field structures in C++:
    DnT_fstruct1D struct and operator<<
    DnT_fstruct2D struct and operator<<
    DnT_fstruct3D struct and operator<<

  Put a fstruct array in a Python list:
    template <typename FS>
    py::list print_fstructarray(py::array_t<FS, 0> arr)

  Copy vector components from a field struct to a double array:
    template <typename FS>
    void fstruct_to_double(FS& fs, double* components)

  Copy vector components from a double array to a field struct:
    template <typename FS>
    void doublef_to_struct(double* components, FS& fs)

  The PYBIND11 interfaces to Python for all the above.

*/

#ifndef FSTRUCT_H
#define FSTRUCT_H

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>

namespace py = pybind11;

// For (x,) particle coordinates
struct DnT_fstruct1D {
// Using PYBIND11_NUMPY_DTYPE_EX, the C++ variable names can be different from the Python names
  double x_;                  // f8    8 bytes
};

// Declare the << operator for the DnT_fstruct1D data type. This is used in print_fstructarray()
std::ostream& operator<<(std::ostream& os, const DnT_fstruct1D& p);

// For (x, y) particle coordinates
struct DnT_fstruct2D {
// Using PYBIND11_NUMPY_DTYPE_EX, the C++ variable names can be different from the Python names
  double x_;                  // f8    8 bytes
  double y_;                  // f8   16
};

// Declare the << operator for the DnT_fstruct2D data type. This is used in print_fstructarray()
std::ostream& operator<<(std::ostream& os, const DnT_fstruct2D& p);

// For (x, y, z) particle coordinates
struct DnT_fstruct3D {
// Using PYBIND11_NUMPY_DTYPE_EX, the C++ variable names can be different from the Python names
  double x_;                  // f8    8 bytes
  double y_;                  // f8   16
  double z_;                  // f8   24
};

// Declare the << operator for the DnT_fstruct3D data type. This is used in print_fstructarray()
std::ostream& operator<<(std::ostream& os, const DnT_fstruct3D& p);

// Declare the C++ function print_fstructarray(): it's templated on the type of the struct
// that corresponds to the Numpy structured array. It takes the array of structs, writes
// all the values into a string, and returns the string. This function is used later to
// define a Python-callable function to write the values in a structured array.
template <typename FS>
py::list print_fstructarray(py::array_t<FS, 0> arr);

// This function copies the spatial coordinates from a field structure array to a double array.
// It uses template specialization to handle field structures with different dimensions.
template <typename FS>
void fstruct_to_double(FS& fs, double* components);

// This function copies a vector of doubles to a field struct.
// It uses template specialization to handle different dimensions.
template <typename FS>
void double_to_fstruct(double* components, FS& fs);

#endif
