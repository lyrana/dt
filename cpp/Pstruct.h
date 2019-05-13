// Copyright (C) 2018 L. D. Hughes

// A Numpy structured array corresponds to a struct in C/C++. The struct
// member-variable names have to be the same as the field names in the Python structured
// array unless PYBIND11_NUMPY_DTYPE_EX is used.  The following struct corresponds to
// a particle structured in Python.

/*
Contents:
  Particle structs in C++:
    dnt::pstruct1D struct and operator<<
    dnt::pstruct2D struct and operator<<
    dnt::pstruct3D struct and operator<<

  Put a pseg array in a Python list:
    template <typename PS>
    py::list print_pstructarray(py::array_t<PS, 0> arr)

  Copy spatial coordinates from a particle struct to a double array:
    template <typename PS>
    void pstruct_to_double(PS& ps, double* point)

  There are currently no PYBIND11 interfaces for the above.

*/

#ifndef PSTRUCT_H
#define PSTRUCT_H

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
// #include <pybind11/eigen.h>
#include <pybind11/operators.h>

namespace py = pybind11;

namespace dnt
{

  /*! \class Ptype
      \brief Ptype is an enum class enumerating the particle structure types.
  */
  enum class Ptype
  {
    cartesian_x,
    cartesian_x_y,
    cartesian_x_y_z
  };
  
/*! \struct Pstruct

  \brief The Pstruct for an unspecified particle structure type.

  \param PT is the particle structure type.
  \sa Ptype, Pstruct_cartesian_x, 

*/
  template <Ptype PT>
  class Pstruct {
    
  public:

    static int DELETE_FLAG;
    static int TRAJECTORY_FLAG;
    
/*    
    push_back(py::tuple ptuple)
    {
      std::cout << "Hello from the set" << std::endl;
    }
*/

/*    
    friend std::ostream& operator<<(std::ostream& os, const Pstruct<PT>& p)
    {
      T ptype;
      
      switch(ptype)
      {
      case Ptype::cartesian_x     : os << "cartesian_x"; break;
      case Ptype::cartesian_x_y   : os << "cartesian_x_y"; break;
      case Ptype::cartesian_x_y_z : os << "cartesian_x_y_z"; break;
      default              : os.setstate(std::ios_base::failbit);
      }
      return os;
    }

*/
  };
  //  class Pstruct ENDCLASS

  // Set the bit patterns for flags. These are static class members, so they're set
  // outside the ctor.
  template <Ptype PT>
    int Pstruct<PT>::DELETE_FLAG = 0b1;  // the lowest bit is 1
  template <Ptype PT>
    int Pstruct<PT>::TRAJECTORY_FLAG = 0b1 << 1; // the second lowest bit is 1
  
/*! \struct Pstruct_cartesian_x

  \brief struct for (x,) particle coordinates.

  \sa Ptype, Pstruct,

*/
  template<>
    class Pstruct<Ptype::cartesian_x>
    {
      // Using PYBIND11_NUMPY_DTYPE_EX, the C++ variable names can be different from the Python names
    public:
      double x_;                  // 1 f8    8 bytes
      double x0_;                 // 2 f8   16 
      double ux_;                 // 3 f8   24
      double weight_;             // 4 f8   32
      int bitflags_;              // 5 i4   36
      int cell_index_;            // 6 i4   40
      int unique_ID_;             // 7 i4   44
      int crossings_;             // 8 i4   48

    public:
      //! Set the member values from a tuple
      void set_from_tuple(py::tuple p)
      {
        x_ = p[0].cast<double>();
        x0_ = p[1].cast<double>();
        ux_ = p[2].cast<double>();
        weight_ = p[3].cast<double>();
        bitflags_ = p[4].cast<int>();
        cell_index_ = p[5].cast<int>();
        unique_ID_ = p[6].cast<int>();
        crossings_ = p[7].cast<int>();
      }

      //! Put the member values into a tuple
      py::tuple as_tuple()
        {
          // Create a tuple containing the member values
          return py::make_tuple(x_, x0_, ux_, weight_, bitflags_, cell_index_, unique_ID_, crossings_);
        }
    
      // Declare the << operator for the pstruct1D data type. This is used in print_pstructarray()
      friend std::ostream& operator<<(std::ostream& os, const Pstruct<Ptype::cartesian_x>& p)
      {
        os << "cartesian_x. fixme";
        return os;
      }
    };
  // class Pstruct<Ptype::cartesian_x> ENDCLASS

/*! \struct pstruct2D
  \brief struct for (x, y) particle coordinates
*/
  template<>
    class Pstruct<Ptype::cartesian_x_y>
    {
      // Using PYBIND11_NUMPY_DTYPE_EX, the C++ variable names can be different from the Python names
    public:
      double x_;                  // 1  f8    8 bytes
      double y_;                  // 2  f8   16
      double x0_;                 // 3  f8   24
      double y0_;                 // 4  f8   32
      double ux_;                 // 5  f8   40
      double uy_;                 // 6  f8   48
      double weight_;             // 7  f8   56
      int bitflags_;              // 8  i4   60
      int cell_index_;            // 9  i4   64
      int unique_ID_;             // 10 i4   68
      int crossings_;             // 11 i4   72

    public:
      //! Set the member values from a tuple
      void set_from_tuple(py::tuple p)
      {
        x_ = p[0].cast<double>();
        y_ = p[1].cast<double>();
        x0_ = p[2].cast<double>();
        y0_ = p[3].cast<double>();
        ux_ = p[4].cast<double>();
        uy_ = p[5].cast<double>();
        weight_ = p[6].cast<double>();
        bitflags_ = p[7].cast<int>();
        cell_index_ = p[8].cast<int>();
        unique_ID_ = p[9].cast<int>();
        crossings_ = p[10].cast<int>();
      }

      //! Put the member values into a tuple
      py::tuple as_tuple()
        {
          // Create a tuple containing the member values
          return py::make_tuple(x_, y_, x0_, y0_, ux_, uy_, weight_, bitflags_, cell_index_, unique_ID_, crossings_);
        }
    
      // Declare the << operator for the pstruct1D data type. This is used in print_pstructarray()
      friend std::ostream& operator<<(std::ostream& os, const Pstruct<Ptype::cartesian_x_y>& p)
      {
        os << "cartesian_x_y. fixme";
        return os;
      }
      
    };
  // class Pstruct<Ptype::cartesian_x_y> ENDCLASS


/*! \struct pstruct3D
  \brief struct for (x, y, z) particle coordinates
*/
  template<>
    class Pstruct<Ptype::cartesian_x_y_z>
  {
// Using PYBIND11_NUMPY_DTYPE_EX, the C++ variable names can be different from the Python names
  public:
    double x_;                  // 1  f8    8 bytes
    double y_;                  // 2  f8   16
    double z_;                  // 3  f8   24
    double x0_;                 // 4  f8   32
    double y0_;                 // 5  f8   40
    double z0_;                 // 6  f8   48
    double ux_;                 // 7  f8   56
    double uy_;                 // 8  f8   64
    double uz_;                 // 9  f8   72
    double weight_;             // 10 f8   80
    int bitflags_;              // 11 i4   84
    int cell_index_;            // 12 i4   88
    int unique_ID_;             // 13 i4   92
    int crossings_;             // 14 i4   96

  public:
      //! Set the member values from a tuple
      void set_from_tuple(py::tuple p)
      {
        x_ = p[0].cast<double>();
        y_ = p[1].cast<double>();
        z_ = p[2].cast<double>();
        x0_ = p[3].cast<double>();
        y0_ = p[4].cast<double>();
        z0_ = p[5].cast<double>();
        ux_ = p[6].cast<double>();
        uy_ = p[7].cast<double>();
        uz_ = p[8].cast<double>();
        weight_ = p[9].cast<double>();
        bitflags_ = p[10].cast<int>();
        cell_index_ = p[11].cast<int>();
        unique_ID_ = p[12].cast<int>();
        crossings_ = p[13].cast<int>();
      }

      //! Put the member values into a tuple
      py::tuple as_tuple()
        {
          // Create a tuple containing the member values
          return py::make_tuple(x_, y_, z_, x0_, y0_, z0_, ux_, uy_, uz_, weight_, bitflags_, cell_index_, unique_ID_, crossings_);
        }
    
      // Declare the << operator for the pstruct1D data type. This is used in print_pstructarray()
      friend std::ostream& operator<<(std::ostream& os, const Pstruct<Ptype::cartesian_x_y_z>& p)
      {
        os << "cartesian_x_y_z. fixme";
        return os;
      }
    
  };
  // class Pstruct<Ptype::cartesian_x_y_z> ENDCLASS

// Declare the << operator for the pstruct3D data type. This is used in print_pstructarray()
//  std::ostream& operator<<(std::ostream& os, const pstruct3D& p);

// Declare the C++ function print_pstructarray(): it's templated on the type of the struct
// that corresponds to the Numpy structured array. It takes the array of structs, writes
// all the values into a string, and returns the string. This function is used later to
// define a Python-callable function to write the values in a structured array.

  
  /* template <typename PS> */
  /*   py::list print_pstructarray(py::array_t<PS, 0> arr); */

// This function copies the spatial coordinates from a particle struct to a double array.
// It uses template specialization to handle particle structs with different dimensions.
  
  /* template <typename PS> */
  /*   void pstruct_to_double(PS& ps, double* point); */

} // namespace dnt

#endif
