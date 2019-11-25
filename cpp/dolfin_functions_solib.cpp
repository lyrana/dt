/*! \file dolfin_functions_solib.cpp

  \brief This file creates a shared library with the Python bindings for
  the C++ dolfin-based functions.

  \namespace dnt
  \sa dolfin_functions.h, dolfin_functions.cpp, segmentedarraypair_solib.cpp

  This file contains the Python-to-C++ bindings allows access to the C++ versions of
  functions defined in dolfin_functions.cpp.

*/

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
  PYBIND11_MODULE(dolfin_functions_solib, m) {

    // C++ functions defined in dolfin_functions.cpp that are called from Python.

    //  C++ signature of is_inside_vertices() is:
    //    bool is_inside_vertices(dolfin::Mesh& mesh, const unsigned int* vertices, double* point);
    m.def("is_inside_vertices", [](dolfin::Mesh& mesh, py::list vertices, py::list point)
          {
            // Create const unsigned int* to hold vertices
            unsigned int verticesTmp[4];
            //1 std::cout << "is_inside_vertices vertices:";
            for (size_t i = 0; i < vertices.size(); i++)
              {
                verticesTmp[i] = vertices[i].cast<unsigned int>();
                //2 std::cout << " " << verticesTmp[i];
              }
            //3 std::cout << std::endl;            
            // Create double* to hold point
            double pointTmp[3];
            for (size_t i = 0; i < point.size(); i++)
              {
                pointTmp[i] = point[i].cast<double>();
              }
            //            return is_inside_vertices(mesh, verticesTmp, pointTmp);
            auto YesNo = is_inside_vertices(mesh, verticesTmp, pointTmp);
            //3 std::cout << "YesNo is " << YesNo << std::endl;
            return YesNo;
          });
    
    //  C++ signature of find_facet() is:
    //    py::tuple find_facet(py::object mesh_M, double* x0, double* dx, size_t cell_index);
    m.def("find_facet", [](py::object mesh_M, py::list x, py::list dx, std::size_t cell_index)
          {
            // Create double* to hold x, dx
            double xTmp[3], dxTmp[3];
            for (size_t i = 0; i < x.size(); i++)
              {
                xTmp[i] = x[i].cast<double>();
                dxTmp[i] = dx[i].cast<double>();
              }
            return find_facet(mesh_M, xTmp, dxTmp, cell_index);
          });
          
    // Call the C++ functions defined in dolfin_functions.cpp
    // This is a check that the functions have been compiled, and have the right
    // argument types. If they haven't, the compiler will give an error.
  
    dolfin::Mesh mesh;
    dolfin::Mesh& meshRef = mesh;
    const unsigned int* vertices(nullptr);
    double* point(nullptr);

    is_inside_vertices(meshRef, vertices, point);

    // Write one for find_facet():
  }

} // namespace dnt
