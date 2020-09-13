// Copyright (C) 2018 L. D. Hughes

// A Numpy structured array corresponds to a array of structs in C/C++. The structure
// member-variable names have to be the same as the field names in the Python structured
// array unless PYBIND11_NUMPY_DTYPE_EX is used.  The following structs corresponds to
// single elements of 1, 2, and 3D structured Numpy arrays of particles in Python.

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
    void pstruct_to_point(PS& ps, double* point)

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
//#include <pybind11/cast.h>

namespace py = pybind11;

namespace dnt
{
  /*! \class Ptype
      \brief Ptype is an enum class enumerating the particle structure types.
  */
  enum class Ptype
  {
    cartesian_x,
    cartesian_xy,
    cartesian_xyz,
    spherical_r
   };

  
/*! \struct Pstruct

  \brief The Pstruct for an unspecified particle structure type.

  \param PT is the particle structure type.
  \sa Ptype, Pstruct_cartesian_x

*/
  template <Ptype PT>
    class Pstruct {

  public:
    // These are static class members, so they're set outside the header file, in
    // particle_solib.cpp
    static int DELETE_FLAG;
    static int TRAJECTORY_FLAG;
    
  };
  //  class Pstruct ENDCLASS


/*! \struct Pstruct_cartesian_x

  \brief struct for (x,) particle coordinates.

  \sa Ptype, Pstruct,

*/
  template<>
    class Pstruct<Ptype::cartesian_x>
    {
      // Using PYBIND11_NUMPY_DTYPE_EX, the C++ variable names can be different from the Python names
    public:
    // These are static class members, so they're set outside the header file, in
    // particle_solib.cpp
      static int DELETE_FLAG;
      static int TRAJECTORY_FLAG;
      
      double x_;                  // 1 f8    8 bytes
      double x0_;                 // 2 f8   16 
      double ux_;                 // 3 f8   24
      double weight_;             // 4 f8   32
      int bitflags_;              // 5 i4   36
      int cell_index_;            // 6 i4   40
      int unique_ID_;             // 7 i4   44
      int crossings_;             // 8 i4   48

    public:

      //! Specularly reflect a particle from a surface.
      /*!

        The surface is specified by the unit normal surface_normal.

        \param[in] surface_normal   Unit normal to the reflecting surface.
        \param[in] dx               Displacement of particle through the surface.

       */
      inline void reflect_from_surface(py::array_t<double> &surface_normal, const double dx[])
      {
        auto sn = surface_normal.unchecked<1>();

        auto n_dot_dx = sn[0] * dx[0];
        auto n_dot_u = sn[0] * ux_;

        x_ -= 2.0 * n_dot_dx * sn[0];
        ux_ -= 2.0 * n_dot_u * sn[0];
      }

      // Overload 2 versions of set_from_list_or_tuple()
      
      //! Set the member values from a py::tuple
      void set_from_list_or_tuple(py::tuple p)
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
      // ENDDEF: void set_from_tuple(py::tuple p)

      //! Set the member values from a py::list
      void set_from_list_or_tuple(py::list p)
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
      // ENDDEF: void set_from_list(py::list p)

      //! Return a particle's data as a tuple
      py::tuple as_tuple()
        {
          // Create a tuple containing the member values
          return py::make_tuple(x_, x0_, ux_, weight_, bitflags_, cell_index_, unique_ID_, crossings_);
        }
      // ENDDEF: py::tuple as_tuple()

      //! Return a particle's data as a list
      py::list as_list()
        {
          // Create a list containing the member values
          auto l = py::list();
          l.append(x_);
          l.append(x0_);
          l.append(ux_);
          l.append(weight_);
          l.append(bitflags_);
          l.append(cell_index_);
          l.append(unique_ID_);
          l.append(crossings_);
          
          return l;
        }
      // ENDDEF: py::list as_list()

      
      //! Return a particle's data as a dictionary
      py::dict as_dict()
        {
          using namespace pybind11::literals;
          // Create a dictionary containing the member values
          return py::dict("x"_a=x_, "x0"_a=x0_, "ux"_a=ux_, "weight"_a=weight_, "bitflags"_a=bitflags_, "cell_index"_a=cell_index_, "unique_ID"_a=unique_ID_, "crossings"_a=crossings_);
        }
      // ENDDEF: py::tuple as_dict()


      // Declare the << operator for the Ptype::cartesian_x data type. This is used in print_pstructarray()
      friend std::ostream& operator<<(std::ostream& os, const Pstruct<Ptype::cartesian_x>& p)
      {
        os << "cartesian_x. fixme";
        return os;
      }
      // ENDDEF: operator<<()
    };
  // ENDCLASS: class Pstruct<Ptype::cartesian_x>

/*! \struct pstruct2D
  \brief struct for (x, y) particle coordinates
*/
  template<>
    class Pstruct<Ptype::cartesian_xy>
    {
      // Using PYBIND11_NUMPY_DTYPE_EX, the C++ variable names can be different from the Python names
    public:
    // These are static class members, so they're set outside the header file, in
    // particle_solib.cpp
      static int DELETE_FLAG;
      static int TRAJECTORY_FLAG;
      
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

      //! Specularly reflect a particle from a surface.
      /*!

        The surface is specified by the unit normal surface_normal.

        \param[in] surface_normal   Unit normal to the reflecting surface.
        \param[in] dx               Displacement of particle through the surface.

       */
      inline void reflect_from_surface(py::array_t<double> &surface_normal, const double dx[])
      {
        auto sn = surface_normal.unchecked<1>();

        auto n_dot_dx = sn[0]*dx[0] + sn[1]*dx[1];
        auto n_dot_u = sn[0]*ux_ + sn[1]*uy_;

        x_ -= 2.0 * n_dot_dx * sn[0];
        y_ -= 2.0 * n_dot_dx * sn[1];
        ux_ -= 2.0 * n_dot_u * sn[0];
        uy_ -= 2.0 * n_dot_u * sn[1];
      }
      
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

      //! Set the member values from a list
      void set_from_list_or_tuple(py::list p)
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
    
      //! Put the member values into a py::list
      py::list as_list()
        {
          // Create a list containing the member values
          auto l = py::list();
          l.append(x_);
          l.append(y_);
          l.append(x0_);
          l.append(y0_);
          l.append(ux_);
          l.append(uy_);
          l.append(weight_);
          l.append(bitflags_);
          l.append(cell_index_);
          l.append(unique_ID_);
          l.append(crossings_);
          
          return l;
        }
      // py::list as_list() ENDDEF
      
      //! Put the member values into a dictionary
      py::dict as_dict()
        {
          using namespace pybind11::literals;          
          // Create a dictionary containing the member values
          return py::dict("x"_a=x_, "y"_a=y_, "x0"_a=x0_, "y0"_a=y0_, "ux"_a=ux_, "uy"_a=uy_, "weight"_a=weight_, "bitflags"_a=bitflags_, "cell_index"_a=cell_index_, "unique_ID"_a=unique_ID_, "crossings"_a=crossings_);
        }
      // py::tuple as_dict() ENDDEF

      // Declare the << operator for the Ptype::cartesian_xy data type. This is used in print_pstructarray()
      friend std::ostream& operator<<(std::ostream& os, const Pstruct<Ptype::cartesian_xy>& p)
      {
        return os <<   "x=" << p.x_
                  << ", y=" << p.y_
                  << ", x0=" << p.x0_
                  << ", y0=" << p.y0_
                  << ", ux=" << p.ux_
                  << ", uy=" << p.uy_
                  << ", weight=" << p.weight_
                  << ", bitflags=" << p.bitflags_
                  << ", cell_index=" << p.cell_index_
                  << ", unique_ID=" << p.unique_ID_
                  << ", crossings=" << p.crossings_;
      }
      
    };
  // class Pstruct<Ptype::cartesian_xy> ENDCLASS


/*! \struct pstruct3D
  \brief struct for (x, y, z) particle coordinates
*/
  template<>
    class Pstruct<Ptype::cartesian_xyz>
  {
// Using PYBIND11_NUMPY_DTYPE_EX, the C++ variable names can be different from the Python names
  public:
    // These are static class members, so they're set outside the header file, in
    // particle_solib.cpp
    static int DELETE_FLAG;
    static int TRAJECTORY_FLAG;
    
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

      //! Specularly reflect a particle from a surface.
      /*!

        The surface is specified by the unit normal surface_normal.

        \param[in] surface_normal   Unit normal to the reflecting surface.
        \param[in] dx               Displacement of particle through the surface.

       */
      inline void reflect_from_surface(py::array_t<double> &surface_normal, const double dx[])
      {
        auto sn = surface_normal.unchecked<1>();

        auto n_dot_dx = sn[0]*dx[0] + sn[1]*dx[1] + sn[2]*dx[2];
        auto n_dot_u = sn[0]*ux_ + sn[1]*uy_ + sn[2]*uz_;

        x_ -= 2.0 * n_dot_dx * sn[0];
        y_ -= 2.0 * n_dot_dx * sn[1];
        z_ -= 2.0 * n_dot_dx * sn[2];
        ux_ -= 2.0 * n_dot_u * sn[0];
        uy_ -= 2.0 * n_dot_u * sn[1];
        uz_ -= 2.0 * n_dot_u * sn[2];
      }
    
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
      //! Set the member values from a list
//      void set_from_list(py::list p)
      void set_from_list_or_tuple(py::list p)
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
    
      //! Put the member values into a py::list
      py::list as_list()
        {
          // Create a list containing the member values
          auto l = py::list();
          l.append(x_);
          l.append(y_);
          l.append(z_);
          l.append(x0_);
          l.append(y0_);
          l.append(z0_);
          l.append(ux_);
          l.append(uy_);
          l.append(uz_);
          l.append(weight_);
          l.append(bitflags_);
          l.append(cell_index_);
          l.append(unique_ID_);
          l.append(crossings_);
          
          return l;
        }
      // py::list as_list() ENDDEF
      
      //! Put the member values into a dictionary
      py::dict as_dict()
        {
          using namespace pybind11::literals;          
          // Create a dictionary containing the member values
          return py::dict("x"_a=x_, "y"_a=y_, "z"_a=z_, "x0"_a=x0_, "y0"_a=y0_, "z0"_a=z0_, "ux"_a=ux_, "uy"_a=uy_, "uz"_a=uz_, "weight"_a=weight_, "bitflags"_a=bitflags_, "cell_index"_a=cell_index_, "unique_ID"_a=unique_ID_, "crossings"_a=crossings_);
        }
      // py::tuple as_dict() ENDDEF

      // Declare the << operator for the Pstruct::cartesian_xyz data type. This is used in print_pstructarray()
      friend std::ostream& operator<<(std::ostream& os, const Pstruct<Ptype::cartesian_xyz>& p)
      {
        os << "Printing a Ptype::cartesian_xyz. fixme";
        return os;
      }
    
  };
  // class Pstruct<Ptype::cartesian_xyz> ENDCLASS

/*! \struct Pstruct_spherical_r

  \brief struct for (r,) particle coordinates.

  \sa Ptype, Pstruct,

*/
  template<>
    class Pstruct<Ptype::spherical_r>
    {
      // Using PYBIND11_NUMPY_DTYPE_EX, the C++ variable names can be different from the Python names
    public:
    // These are static class members, so they're set outside the header file, in
    // particle_solib.cpp
      static int DELETE_FLAG;
      static int TRAJECTORY_FLAG;
      
      double r_;                  // 1 f8    8 bytes
      double r0_;                 // 2 f8   16 
      double ur_;                 // 3 f8   24
      double weight_;             // 4 f8   32
      int bitflags_;              // 5 i4   36
      int cell_index_;            // 6 i4   40
      int unique_ID_;             // 7 i4   44
      int crossings_;             // 8 i4   48

    public:

      //! Specularly reflect a particle from a surface.
      /*!

        The surface is specified by the unit normal surface_normal.

        \param[in] surface_normal   Unit normal to the reflecting surface.
        \param[in] dx               Displacement of particle through the surface.

       */
      inline void reflect_from_surface(py::array_t<double> &surface_normal, const double dx[])
      {
        auto sn = surface_normal.unchecked<1>();

        auto n_dot_dx = sn[0] * dx[0];
        auto n_dot_u = sn[0] * ur_;

        r_ -= 2.0 * n_dot_dx * sn[0];
        ur_ -= 2.0 * n_dot_u * sn[0];
      }

      // Overload 2 versions of set_from_list_or_tuple()
      
      //! Set the member values from a py::tuple
      void set_from_list_or_tuple(py::tuple p)
      {
        r_ = p[0].cast<double>();
        r0_ = p[1].cast<double>();
        ur_ = p[2].cast<double>();
        weight_ = p[3].cast<double>();
        bitflags_ = p[4].cast<int>();
        cell_index_ = p[5].cast<int>();
        unique_ID_ = p[6].cast<int>();
        crossings_ = p[7].cast<int>();
      }
      // ENDDEF: void set_from_tuple(py::tuple p)

      //! Set the member values from a py::list
      void set_from_list_or_tuple(py::list p)
      {
        r_ = p[0].cast<double>();
        r0_ = p[1].cast<double>();
        ur_ = p[2].cast<double>();
        weight_ = p[3].cast<double>();
        bitflags_ = p[4].cast<int>();
        cell_index_ = p[5].cast<int>();
        unique_ID_ = p[6].cast<int>();
        crossings_ = p[7].cast<int>();
      }
      // ENDDEF: void set_from_list(py::list p)

      //! Return a particle's data as a tuple
      py::tuple as_tuple()
        {
          // Create a tuple containing the member values
          return py::make_tuple(r_, r0_, ur_, weight_, bitflags_, cell_index_, unique_ID_, crossings_);
        }
      // ENDDEF: py::tuple as_tuple()

      //! Return a particle's data as a list
      py::list as_list()
        {
          // Create a list containing the member values
          auto l = py::list();
          l.append(r_);
          l.append(r0_);
          l.append(ur_);
          l.append(weight_);
          l.append(bitflags_);
          l.append(cell_index_);
          l.append(unique_ID_);
          l.append(crossings_);
          
          return l;
        }
      // ENDDEF: py::list as_list()

      
      //! Return a particle's data as a dictionary
      py::dict as_dict()
        {
          using namespace pybind11::literals;
          // Create a dictionary containing the member values
          return py::dict("r"_a=r_, "r0"_a=r0_, "ur"_a=ur_, "weight"_a=weight_, "bitflags"_a=bitflags_, "cell_index"_a=cell_index_, "unique_ID"_a=unique_ID_, "crossings"_a=crossings_);
        }
      // ENDDEF: py::tuple as_dict()


      // Declare the << operator for the Ptype::spherical_r data type.
      friend std::ostream& operator<<(std::ostream& os, const Pstruct<Ptype::spherical_r>& p)
      {
        return os <<   "r=" << p.r_
                  << ", r0=" << p.r0_
                  << ", ur=" << p.ur_
                  << ", weight=" << p.weight_
                  << ", bitflags=" << p.bitflags_
                  << ", cell_index=" << p.cell_index_
                  << ", unique_ID=" << p.unique_ID_
                  << ", crossings=" << p.crossings_;
      }
      // ENDDEF: operator<<()

    };
  // ENDCLASS: class Pstruct<Ptype::spherical_r>

  
// Declare the C++ function print_pstructarray(): it's templated on the type of the struct
// that corresponds to the Numpy structured array. It takes the array of structs, writes
// all the values into a string, and returns the string. This function is used later to
// define a Python-callable function to write the values in a structured array.
  
  /* template <typename PS> */
  /*   py::list print_pstructarray(py::array_t<PS, 0> arr); */

// This function copies the spatial coordinates from a particle struct to a double array.
// It uses template specialization to handle particle structs with different dimensions.
  
  /* template <typename PS> */
  /*   void pstruct_to_point(PS& ps, double* point); */

} // namespace dnt

#endif
