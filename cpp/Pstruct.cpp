
#include "pstruct.h"

namespace py = pybind11;

/*
void init_ex1(py::module &m) {
m.def("add", [](int a, int b) { return a + b; });
}
*/

namespace dnt
{
// Define the << operator for the DnT_pstruct1D data type. This is used in print_pstructarray()
// below to print each element in the array.
  std::ostream& operator<<(std::ostream& os, const Pstruct<Ptype::cartesian_x>& p)
  //  std::ostream& operator<<(std::ostream& os, const Ptype& p)
  {
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
  std::ostream& operator<<(std::ostream& os, const Pstruct<Ptype::cartesian_x_y>& p) {
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
  std::ostream& operator<<(std::ostream& os, const Pstruct<Ptype::cartesian_x_y_z>& p) {
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
  template <Ptype PT>
  py::list print_pstructarray(py::array_t<Pstruct<PT>, 0> arr) {
    const auto arr_info = arr.request();  // request() returns metadata about the array (ptr, ndim, size, shape)
    const auto p = static_cast<Pstruct<PT>*>(arr_info.ptr); // Pointer to a particle structure in pseg

// See ~/local/include/pybind11/buffer_info.h
// arr_info attributes are:
// itemsize, size, format, ndim, shape, strides 
    std::cout << "pstruct.cpp:print_pstructarray():"
              << " itemsize=" << arr_info.itemsize
              << ", size=" << arr_info.size
              << ", format=" << arr_info.format
              << ", ndim=" << arr_info.ndim
              << std::endl;

    std::cout << "pstruct.cpp:print_pstructarray(): The size one structure is " << sizeof(Pstruct<PT>) << " bytes" << std::endl;
  
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
  template <Ptype PT>
  void pstruct_to_double(Pstruct<PT>& p, double* point)
  {
    // There is no generic version
  }
  template <>
  void pstruct_to_double<Pstruct<Ptype::cartesian_x>(Pstruct<Ptype::cartesian_x& p, double* point)
  {
    point[0] = p.x_;
    // how about point = &(p.x_);
   
  }
  template <>
  void pstruct_to_double<Pstruct<Ptype::cartesian_x_y>(Pstruct<Ptype::cartesian_x_y& p, double* point)
  {
    point[0] = p.x_;
    point[1] = p.y_;
  }
  template <>
  void pstruct_to_double<Pstruct<Ptype::cartesian_x_y_z>(Pstruct<Ptype::cartesian_x_y_z& p, double* point)
  {
    point[0] = p.x_;
    point[1] = p.y_;
    point[2] = p.z_;
  }

//! Force the compiler to generate the specialized forms of the templated function print_pstructarray()

/*!
  The templated signature is:
  template <typename S>
  py::list print_pstructarray(py::array_t<Pstruct<PT>, 0> arr) {
*/

  void force_DnT_print_pstructarray_instances()
  {
    py::array_t<DnT_pstruct1D, 0> pseg1D;
    py::array_t<DnT_pstruct2D, 0> pseg2D;
    py::array_t<DnT_pstruct3D, 0> pseg3D;

    print_pstructarray<DnT_pstruct1D>(pseg1D);
    print_pstructarray<DnT_pstruct2D>(pseg2D);
    print_pstructarray<DnT_pstruct3D>(pseg3D);
  
  }
}
