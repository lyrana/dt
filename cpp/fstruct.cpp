
#include "fstruct.h"

namespace py = pybind11;

/*
void init_ex1(py::module &m) {
m.def("add", [](int a, int b) { return a + b; });
}
*/


// Define the << operator for the DnT_fstruct1D data type. This is used in print_fstructarray()
// below to print each element in the array.
std::ostream& operator<<(std::ostream& os, const DnT_fstruct1D& fs) {
  return os <<   "fx=" << fs.x_;
};

// Define the << operator for the DnT_fstruct2D data type. This is used in print_fstructarray()
// below to print each element in the array.
std::ostream& operator<<(std::ostream& os, const DnT_fstruct2D& fs) {
  return os <<   "fx=" << fs.x_
            << ", fy=" << fs.y_;
};

// Define the << operator for the DnT_fstruct3D data type. This is used in print_fstructarray()
// below to print each element in the array.
std::ostream& operator<<(std::ostream& os, const DnT_fstruct3D& fs) {
  return os <<   "fx=" << fs.x_
            << ", fy=" << fs.y_
            << ", fz=" << fs.z_;
};



// Define the C++ function print_fstructarray(): it's templated on the type of the struct
// that corresponds to the Numpy structured array. It takes the array of structs, writes
// all the values into a string, and returns the string. This function is used later to
// define a Python-callable function to write the values in a structured array.
template <typename FS>
py::list print_fstructarray(py::array_t<FS, 0> arr) {
  const auto arr_info = arr.request();  // request() returns metadata about the array (ptr, ndim, size, shape)
  const auto fs = static_cast<FS*>(arr_info.ptr); // Pointer to a field struct in pseg

// See ~/local/include/pybind11/buffer_info.h
// arr_info attributes are:
// itemsize, size, format, ndim, shape, strides 
  std::cout << "fstruct.cpp:print_fstructarray(): \n"
            << "  itemsize=" << arr_info.itemsize
            << ", size=" << arr_info.size
            << ", format=" << arr_info.format
            << ", ndim=" << arr_info.ndim
            << std::endl;

  std::cout << "fstruct.cpp:print_fstructarray(): \n The size one array element is " << sizeof(FS) << " bytes" << std::endl;
  
  auto l = py::list();
  for (ssize_t i = 0; i < arr_info.size; i++) {
    std::stringstream ss;
    ss << fs[i];  // Use the above definition of "<<" to write the i'th struct as a string.
    l.append(py::str(ss.str()));
  }
  return l;
}

// This function copies values from a field struct to a vector of doubles.
// It uses template specialization to handle different dimensions.
template <typename FS>
void fstruct_to_double(FS& fs, double* components)
{
  // There is no generic version
}
template <>
void fstruct_to_double<DnT_fstruct1D>(DnT_fstruct1D& fs1D, double* components)
{
  components[0] = fs1D.x_;
  // how about components = &(fs1D.x_); ?
   
}
template <>
void fstruct_to_double<DnT_fstruct2D>(DnT_fstruct2D& fs2D, double* components)
{
  components[0] = fs2D.x_;
  components[1] = fs2D.y_;
}
template <>
void fstruct_to_double<DnT_fstruct3D>(DnT_fstruct3D& fs3D, double* components)
{
  components[0] = fs3D.x_;
  components[1] = fs3D.y_;
  components[2] = fs3D.z_;
}

// This function copies a vector of doubles to a field struct.
// It uses template specialization to handle different dimensions.
template <typename FS>
void double_to_fstruct(double* components, FS& fs)
{
  // There is no generic version
}
template <>
void double_to_fstruct<DnT_fstruct1D>(double* components, DnT_fstruct1D& fs1D)
{
  fs1D.x_ = components[0];
  // how about components = &(fs1D.x_); ?
   
}
template <>
void double_to_fstruct<DnT_fstruct2D>(double* components, DnT_fstruct2D& fs2D)
{
  fs2D.x_ = components[0];
  fs2D.y_ = components[1];
}
template <>
void double_to_fstruct<DnT_fstruct3D>(double* components, DnT_fstruct3D& fs3D)
{
  fs3D.x_ = components[0];
  fs3D.y_ = components[1];
  fs3D.z_ = components[2];
}
