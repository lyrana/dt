/*! \file dolfin_functions_solib.cpp

  \brief This file creates a shared library with the Python bindings for
  the C++ dolfin-based functions in dolfin_functions.cpp.

  These functions are
      find_facet()

      These have two versions: (i) For a Numpy array of points and (ii) for a Pstruct<PT>* with points:
        interpolate_field_to_points_cartesian_xyz()
        interpolate_field_to_points_cartesian_xy()
        interpolate_field_to_points_cartesian_x()

      is_inside_vertices()

  \namespace dnt
  \sa dolfin_functions.h, dolfin_functions.cpp, segmented_array_pair_solib.cpp

  This file contains the Python-to-C++ bindings allows access to the C++ versions of
  functions defined in dolfin_functions.cpp.

*/

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#include "dolfin_functions.h"

namespace py = pybind11;

// The PYBIND11_MODULE() macro creates a function that will be called when an import
// statement is issued from within Python. The module name (dolfin_functions_solib)
// is given as the first macro argument (it should not be in quotes). The second
// argument (m) defines a variable of type py::module which is the main interface for
// creating bindings. The method module::def() generates binding code that exposes
// the C++ functions to Python.

namespace dnt {

  // Create a variable 'm' of type py::module
  // MODULE_NAME can be specified using -DMODULE_NAME= in the makefile.
  PYBIND11_MODULE(MODULE_NAME, m) {

    // m.def makes the Python binding for each C++ function defined in
    // dolfin_functions.cpp.

    m.def("find_facet", [](py::object mesh_M, py::list x0, py::list dx, std::size_t cell_index, bool returnStdArray = false)
          {
            // Convert the above Python args to the args of the C++ function.
            //  C++ signature of find_facet() is:
            //    py::tuple find_facet(py::object mesh_M, double* x0, double* dx, size_t cell_index, bool returnStdArray)
            // Create double* to hold x0, dx. 
            double xTmp[3], dxTmp[3];
            for (size_t i = 0; i < x0.size(); i++)
              {
                xTmp[i] = x0[i].cast<double>();
                dxTmp[i] = dx[i].cast<double>();
              }
            return find_facet(mesh_M, xTmp, dxTmp, cell_index, returnStdArray);
          }, py::arg("mesh_M"), py::arg("x0"), py::arg("dx"), py::arg("cell_index"), py::arg("returnStdArray") = false); // Sets a default value for the last arg.

    // Don't need the PARTICLE_TYPE in the function name if it's in the module name?
    /*
    m.def("interpolate_field_to_points_" TOSTRING(PARTICLE_TYPE), [](py::object field_F, py::array_t<Pstruct<Ptype::PARTICLE_TYPE>, 0> points, py::ssize_t npoints, py::array_t<double> field_at_points)
          {
            auto fieldFunction = field_F.attr("function").cast<dolfin::Function&>();
            interpolate_field_to_points<Ptype::PARTICLE_TYPE>(fieldFunction, points, npoints, field_at_points);
          });
    */

    //
    // Versions for py::array_t<Pstruct<Ptype::PARTICLE_TYPE>, 0> points, i.e., points is a Numpy array.
    //
    m.def("interpolate_field_to_points_cartesian_xyz", [](py::object field_F, py::array_t<Pstruct<Ptype::cartesian_xyz>, 0> points, py::ssize_t npoints, py::array_t<double> field_at_points)
          {
            auto fieldFunction = field_F.attr("function").attr("_cpp_object").cast<dolfin::Function*>();
            // auto fieldFunction = field_F.cast<std::shared_ptr<dolfin::GenericVector>>();
            // std::cout << "fieldFunction->size() = " << fieldFunction->size() << std::endl;
            interpolate_field_to_points<Ptype::cartesian_xyz>(fieldFunction, points, npoints, field_at_points);
            
          });

    m.def("interpolate_field_to_points_cartesian_xy", [](py::object field_F, py::array_t<Pstruct<Ptype::cartesian_xy>, 0> points, py::ssize_t npoints, py::array_t<double> field_at_points)
          {
            auto fieldFunction = field_F.attr("function").attr("_cpp_object").cast<dolfin::Function*>();
            // auto fieldFunction = field_F.cast<std::shared_ptr<dolfin::GenericVector>>();
            // std::cout << "fieldFunction->size() = " << fieldFunction->size() << std::endl;
            interpolate_field_to_points<Ptype::cartesian_xy>(fieldFunction, points, npoints, field_at_points);
            
          });

    m.def("interpolate_field_to_points_cartesian_x", [](py::object field_F, py::array_t<Pstruct<Ptype::cartesian_x>, 0> points, py::ssize_t npoints, py::array_t<double> field_at_points)
          {
            auto fieldFunction = field_F.attr("function").attr("_cpp_object").cast<dolfin::Function*>();
            // auto fieldFunction = field_F.cast<std::shared_ptr<dolfin::GenericVector>>();
            // std::cout << "fieldFunction->size() = " << fieldFunction->size() << std::endl;
            interpolate_field_to_points<Ptype::cartesian_x>(fieldFunction, points, npoints, field_at_points);

          });
    //
    // Versions for Pstruct<PT>* points, i.e., points is a pointer to the particle-structs array.
    //
    m.def("interpolate_field_to_points_cartesian_xyz", [](py::object field_F, Pstruct<Ptype::cartesian_xyz>* points, py::ssize_t npoints, py::array_t<double> field_at_points)
          {
            auto fieldFunction = field_F.attr("function").attr("_cpp_object").cast<dolfin::Function*>();
            // auto fieldFunction = field_F.cast<std::shared_ptr<dolfin::GenericVector>>();
            // std::cout << "fieldFunction->size() = " << fieldFunction->size() << std::endl;
            interpolate_field_to_points<Ptype::cartesian_xyz>(fieldFunction, points, npoints, field_at_points);
            
          });
    m.def("interpolate_field_to_points_cartesian_xy", [](py::object field_F, Pstruct<Ptype::cartesian_xy>* points, py::ssize_t npoints, py::array_t<double> field_at_points)
          {
            auto fieldFunction = field_F.attr("function").attr("_cpp_object").cast<dolfin::Function*>();
            // auto fieldFunction = field_F.cast<std::shared_ptr<dolfin::GenericVector>>();
            // std::cout << "fieldFunction->size() = " << fieldFunction->size() << std::endl;
            interpolate_field_to_points<Ptype::cartesian_xy>(fieldFunction, points, npoints, field_at_points);
            
          });

    m.def("interpolate_field_to_points_cartesian_x", [](py::object field_F, Pstruct<Ptype::cartesian_x>* points, py::ssize_t npoints, py::array_t<double> field_at_points)
          {
            auto fieldFunction = field_F.attr("function").attr("_cpp_object").cast<dolfin::Function*>();
            // auto fieldFunction = field_F.cast<std::shared_ptr<dolfin::GenericVector>>();
            // std::cout << "fieldFunction->size() = " << fieldFunction->size() << std::endl;
            interpolate_field_to_points<Ptype::cartesian_x>(fieldFunction, points, npoints, field_at_points);
            
          });
    
    m.def("is_inside_vertices", [](dolfin::Mesh& mesh, py::list vertices, py::list point)
          {
            // Convert the above Python args to the args of the C++ function
            //  C++ signature of is_inside_vertices() is:
            //    bool is_inside_vertices(dolfin::Mesh& mesh, const unsigned int* vertices, double* point);
            
            // Create unsigned int* to hold vertices
            unsigned int verticesTmp[4];
            for (size_t i = 0; i < vertices.size(); i++)
              {
                verticesTmp[i] = vertices[i].cast<unsigned int>();
              }
            // Create a double* to hold point
            double pointTmp[3];
            for (size_t i = 0; i < point.size(); i++)
              {
                pointTmp[i] = point[i].cast<double>();
              }
            //            return is_inside_vertices(mesh, verticesTmp, pointTmp);
            auto YesNo = is_inside_vertices(mesh, verticesTmp, pointTmp);
            return YesNo;
          });
    
          
    // Call the C++ functions defined in dolfin_functions.cpp
    // This is a check that the functions have been compiled, and have the right
    // argument types. If they haven't, the compiler will give an error.

    // Not clear that these serve any useful purpose. They don't cause a compile of the functions.
    // They just check if the args match the declarations in "dolfin_functions.h" above.

    /*
    dolfin::Mesh mesh;
    dolfin::Mesh& meshRef = mesh;
    const unsigned int* vertices(nullptr);
    double* point(nullptr);

    is_inside_vertices(meshRef, vertices, point);

    dolfin::Function field;
    dolfin::Function& fieldRef = field;
    py::array_t<Pstruct<Ptype::cartesian_xyz>, 0> points;
    py::ssize_t npoints = 0;
    py::array_t<double> field_at_points;

    // This runs when the .so is imported, and bombs:
    interpolate_field_to_points(fieldRef, points, npoints, field_at_points);
    */
    
    // Write one for find_facet():
  }

} // namespace dnt
