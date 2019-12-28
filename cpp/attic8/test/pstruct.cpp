
#include "pstruct.h"

namespace py = pybind11;

/*
void init_ex1(py::module &m) {
m.def("add", [](int a, int b) { return a + b; });
}
*/


// Define the << operator for the DnT_pstruct1D data type. This is used in print_pstructarray()
// below to print each element in the array.
std::ostream& operator<<(std::ostream& os, const DnT_pstruct1D& p) {
  return os <<   "x=" << p.x_
            << ", x0=" << p.x0_
            << ", ux=" << p.ux_
            << ", w=" << p.weight_
            << ", flgs=" << p.bitflags_
            << ", indx=" << p.cell_index_
            << ", id=" << p.unique_ID_
            << ", crossings=" << p.crossings_;
};

// Define the << operator for the DnT_pstruct2D data type. This is used in print_pstructarray()
// below to print each element in the array.
std::ostream& operator<<(std::ostream& os, const DnT_pstruct2D& p) {
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

// Define the << operator for the DnT_pstruct3D data type. This is used in print_pstructarray()
// below to print each element in the array.
std::ostream& operator<<(std::ostream& os, const DnT_pstruct3D& p) {
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



// Define the C++ function print_pstructarray(): it's templated on the type of the struct
// that corresponds to the Numpy structured array. It takes the array of structs, writes
// all the values into a string, and returns the string. This function is used later to
// define a Python-callable function to write the values in a structured array.
template <typename S>
py::list print_pstructarray(py::array_t<S, 0> arr) {
  const auto arr_info = arr.request();  // request() returns metadata about the array (ptr, ndim, size, shape)
  const auto p = static_cast<S*>(arr_info.ptr); // Pointer to a particle structure in pseg

// See ~/local/include/pybind11/buffer_info.h
// arr_info attributes are:
// itemsize, size, format, ndim, shape, strides 
  std::cout << "  itemsize=" << arr_info.itemsize
            << ", size=" << arr_info.size
            << ", format=" << arr_info.format
            << ", ndim=" << arr_info.ndim
            << std::endl;

  std::cout << "The size one structure is " << sizeof(S) << " bytes" << std::endl;
  
  auto l = py::list();
  for (ssize_t i = 0; i < arr_info.size; i++) {
    std::stringstream ss;
    ss << p[i];  // Use the above definition of "<<" to write the i'th structure as a string.
    l.append(py::str(ss.str()));
  }
  return l;
}

// This function copies the spatial coordinates from a particle structure to a double array.
// It uses template specialization to handle particle structures with different dimensions.
template <typename S>
void pstruct_to_double(S& p3D, double* point)
{
  // There is no generic version
}
template <>
void pstruct_to_double<DnT_pstruct1D>(DnT_pstruct1D& p1D, double* point)
{
  point[0] = p1D.x_;
  // how about point = &(p1D.x_);
   
}
template <>
void pstruct_to_double<DnT_pstruct2D>(DnT_pstruct2D& p2D, double* point)
{
  point[0] = p2D.x_;
  point[1] = p2D.y_;
}
template <>
void pstruct_to_double<DnT_pstruct3D>(DnT_pstruct3D& p3D, double* point)
{
  point[0] = p3D.x_;
  point[1] = p3D.y_;
  point[2] = p3D.z_;
}

