#include "dolfin.h"
#include "predicates.h"

namespace py = pybind11;

/*
void init_ex2(py::module &m) {
m.def("sub", [](int a, int b) { return a - b; });
}
*/

namespace dnt
{

  //! Divide the values of a cell function by the cell volumes.
  /*!
    \param dF is a dolfin function defined on the cells of a mesh. dF.vector() is a
    std::shared_ptr<GenericVector>.

    \param cell_volume_dict is a reference to a python dict containing the cell volumes.

    \return void
  */
  void divide_by_cell_volumes(dolfin::Function& dF, const std::map<int, double> &cell_volume_dict)
  {
    std::shared_ptr< const dolfin::FunctionSpace > functionSpace = dF.function_space();
    std::size_t dofs_per_cell = functionSpace->element()->space_dimension();

    // In the case of the cell density, there's only one DoF in the cell, so we use
    // the name 'value' instead of 'values':

    std::vector<double> value(dofs_per_cell);
  
    std::shared_ptr< const dolfin::GenericDofMap > dofMap = functionSpace->dofmap();
      
    auto mesh = functionSpace->mesh();
    // Loop on cells
    /// A CellIterator is a MeshEntityIterator of topological codimension 0.
    /// typedef MeshEntityIteratorBase<Cell> CellIterator;
    for (dolfin::CellIterator cell(*mesh); !cell.end(); ++cell)
      {
        // Get the dofIndex of the cell
        auto cellIndex = cell->index();
        // If there's only one DoF in the cell, so we use the
        // name 'dofIndex' instead of 'dofIndices':
        auto dofIndex = dofMap->cell_dofs(cellIndex);
        dF.vector()->get_local(value.data(), dofs_per_cell, dofIndex.data());
      
        // Get the volume of the cell
        auto cellVol = cell_volume_dict.at(cellIndex);

        // This gives type "d" (double?)
        // std::cout << "type of cellVol: " << typeid(cellVol).name() << std::endl;
      
        // Divide the value by the volume
        value[0] /= cellVol;
        // Put back the value
        dF.vector()->set_local(value.data(), dofs_per_cell, dofIndex.data());        
      }
  }
  // ENDDEF: void divide_by_cell_volumes(dolfin::Function& dF, const std::map<int, double> &cell_volume_dict)  

  //! Evaluate a dolfin Function at a set of points. The point data are in a Python structured array.
  /*!

    This is the C++ version of the Python function in Dolfin_Module.py:
    interpolate_field_to_points(self, points, field_at_points)

    \param field is a reference to a Dolfin Function that has the values of the
    field on a finite-element mesh.

    \param points is a structured Numpy array containing the point data.

    \param field_at_points is a structured Numpy array to hold the values computed at the
    points.

  */
  template <Ptype PT, typename FS>
  void interpolate_field_to_points(dolfin::Function& field,
                                   py::array_t<Pstruct<PT>, 0> points,
                                   py::array_t<FS, 0> field_at_points)
  {
    const auto points_info = points.request(); // request() returns metadata about the Python array (ptr, ndim, size, shape)
    const auto p = static_cast<Pstruct<PT>*>(points_info.ptr); // Pointer to a struct in points
  
    std::shared_ptr< const dolfin::FunctionSpace > functionSpace = field.function_space();
    // This should be 1, since there's just one E value per cell:
    // std::size_t dofs_per_cell = functionSpace->element()->space_dimension();
    std::shared_ptr< const dolfin::GenericDofMap > dofMap = functionSpace->dofmap();
  
    // Temporary variables to access cell information
    const dolfin::Mesh& mesh = *functionSpace->mesh();

    double values[] = {0.0, 0.0, 0.0}; // This holds the returned function values. Could the values be put directly into field_at_points[]?
    for (auto ip = 0; ip < points_info.size; ip++) {
      // std::cout << ip << " x=" << p[ip].x_ << std::endl;

      auto cellIndex = p[ip].cell_index_;
      // std::cout << ip << " cellIndex=" << cellIndex << std::endl;

      // In the case of constant-E-in-cell, there's only one DoF in the cell, so we use the
      // name 'dofIndex' instead of 'dofIndices':
      auto dofIndex = dofMap->cell_dofs(cellIndex);

      // Pull out the values of the (constant) field in this cell
      field.vector()->get_local(values, dofIndex.size(), dofIndex.data());

      // Copy the field vector values to the field struct
      double_to_fstruct<FS>(values, field_at_points[ip]);
    
    }

  }
  //ENDDEF: void interpolate_field_to_points(dolfin::Function& field,...  

  template <Ptype PT, typename FS>
  void interpolate_field_to_points_v1(dolfin::Function& field,
                                      py::array_t<Pstruct<PT>, 0> points,
                                      py::array_t<FS, 0> field_at_points)
  {
    const auto points_info = points.request(); // request() returns metadata about the Python array (ptr, ndim, size, shape)
    const auto p = static_cast<Pstruct<PT>*>(points_info.ptr); // Pointer to a struct in points
  
    std::shared_ptr< const dolfin::FunctionSpace > functionSpace = field.function_space();
  
    // Temporary variables to access cell information
    const dolfin::Mesh& mesh = *functionSpace->mesh();
    ufc::cell ufc_cell;

    double point[] = {0.0, 0.0, 0.0}; // This holds a temporary copy of the point coordinates. Can the point coords be accessed directly?
    double values[] = {0.0, 0.0, 0.0}; // This holds the returned function values. Could the values be put directly into field_at_points[]?
    for (auto ip = 0; ip < points_info.size; ip++) {
      // std::cout << ip << " x=" << p[ip].x_ << std::endl;

      // Copy spatial coordinates of p[ip] to a double[] for passing to evaluate_basis()
      // ?? Can p[ip] be passed directly?
      // Or write a special form of eval!
      pstruct_to_point<Pstruct<PT>>(p[ip], point);
    
      auto cellIndex = p[ip].cell_index_;
      // std::cout << ip << " cellIndex=" << cellIndex << std::endl;

      // Get the cell object containing this particle
      dolfin::Cell cell(mesh, static_cast<std::size_t>(cellIndex));
      cell.get_cell_data(ufc_cell);

      // Evaluate the function at the point
      // See dolfin/dolfin/function/Function.cpp


      // This version of eval() uses the Eigen signature:
    
      //    void Function::eval(Eigen::Ref<Eigen::VectorXd> values,
      //                        Eigen::Ref<const Eigen::VectorXd> x,
      //                        const Cell& dolfin_cell, const ufc::cell& ufc_cell) const
    
      // why does this not work?    
      //    field.eval(Eigen::Map<Eigen::VectorXd> mvalues(values, 3), Eigen::Map<const Eigen::VectorXd> mpoint(point, 3), cell, ufc_cell);
      // Because the "Eigen::Map<Eigen::VectorXd> mvalues(values, 3)" has to be constructed in a separate statement (see following).

      Eigen::Map<Eigen::VectorXd> mvalues(values, 3);    
      Eigen::Map<const Eigen::VectorXd> mpoint(point, 3);
      field.eval(mvalues, mpoint, cell, ufc_cell);

      // This version of eval() uses the dolfin::Array signature
    
      //    void Function::eval(Array<double>& values, const Array<double>& x,
      //                    const Cell& dolfin_cell, const ufc::cell& ufc_cell) const

      // Copy the values to field_at_points (which is a py::array_t<double>)
      // The vector (a rank-1 tensor) has field.value_dimension() components

      //TODO: Can field_at_point[] be used directly, instead of copying?

      double_to_fstruct<FS>(values, field_at_points[ip]);
    
      /*    
            for (std::size_t i = 0; i < field.value_dimension(0); ++i)
            {
            field_at_points[ip];
            }
      */

    }

  }
  // ENDDEF: void interpolate_field_to_points_v1(dolfin::Function& field,...

  //! Test whether a point lies within the cell defined by the vertices.
  /*!
    Overload on the number of particle coordinates: 1D, 2D, 3D 
  */

  /// Versions using Pstruct<>

  // 1D version for Ptype::cartesian_x
  bool cell_contains_point(dolfin::Mesh& mesh,
                           py::array_t<int> vertices,
                           Pstruct<Ptype::cartesian_x> ps1D)
  {
    auto v = vertices.unchecked<1>(); // <1> is the number of array indices. See pybind11 docs.
  
    auto v0 = v(0);
    auto v1 = v(1);

    //  std::cout << "v0: " << v0 << " v1: " << v1 << std::endl;
  
    double p0 = mesh.coordinates()[v0];
    double p1 = mesh.coordinates()[v1];
  
    if (p0 > p1)
      std::swap(p0, p1);
    return p0 <= ps1D.x_ and ps1D.x_ <= p1;
  }
  // ENDDEF: bool cell_contains_point(dolfin::Mesh& mesh,...

  // 2D version for Ptype::cartesian_xy
  // Based on dolfin/dolfin/geometry/CollisionPredicates.cpp: collides_triangle_point_2d()

  bool cell_contains_point(dolfin::Mesh& mesh,
                           py::array_t<int> vertices,
                           Pstruct<Ptype::cartesian_xy> ps2D)
  {
    auto v = vertices.unchecked<1>(); // <1> is the number of array indices. See pybind11 docs.
  
    auto v0 = v(0);
    auto v1 = v(1);
    auto v2 = v(2);

    //  std::cout << "v0: " << v0 << " v1: " << v1 << std::endl;
  
    int meshDimension = 2;
    // mesh.coordinates() returns a 1D array of doubles: x0, y0, x1, y1,... in this
    // case.  So the value of x for vertex n is at n*2.
    double* p0 = &mesh.coordinates()[v0*meshDimension];
    double* p1 = &mesh.coordinates()[v1*meshDimension];
    double* p2 = &mesh.coordinates()[v2*meshDimension];

    // Copy the struct coordinates to a double array
    double point[] = {0.0, 0.0};
    pstruct_to_point<Ptype::cartesian_xy>(ps2D, point);

    // The rest is code copied from CollisionPredicates.cpp
  
    const double ref = dnt::_orient2d(p0, p1, p2);

    if (ref > 0.0)
      {
        return (dnt::_orient2d(p1, p2, point) >= 0.0 and
                dnt::_orient2d(p2, p0, point) >= 0.0 and
                dnt::_orient2d(p0, p1, point) >= 0.0);
      }
    else if (ref < 0.0)
      {
        return (dnt::_orient2d(p1, p2, point) <= 0.0 and
                dnt::_orient2d(p2, p0, point) <= 0.0 and
                dnt::_orient2d(p0, p1, point) <= 0.0);
      }
    else
      {
        return ((dnt::_orient2d(p0, p1, point) == 0.0 and
                 dnt::_collides_segment_point_1d(p0[0], p1[0], point[0]) and
                 dnt::_collides_segment_point_1d(p0[1], p1[1], point[1])) or
                (dnt::_orient2d(p1, p2, point) == 0.0 and
                 dnt::_collides_segment_point_1d(p1[0], p2[0], point[0]) and
                 dnt::_collides_segment_point_1d(p1[1], p2[1], point[1])) or
                (dnt::_orient2d(p2, p0, point) == 0.0 and
                 dnt::_collides_segment_point_1d(p2[0], p0[0], point[0]) and
                 dnt::_collides_segment_point_1d(p2[1], p0[1], point[1])));
      }

  }
  // ENDDEF: bool cell_contains_point(dolfin::Mesh& mesh,...

  // NEW VERSION
  // The call is
  //            cell_contains_point(mesh, vertices,  psegOut[ipOut]);

  bool cell_contains_point(dolfin::Mesh& mesh,
                           py::array_t<int> vertices,
                           double* point)
  {
    // The algorithm used depends on the mesh dimension
    const std::size_t tdim = mesh.topology().dim();
    const std::size_t gdim = mesh.geometry().dim();

    auto v = vertices.unchecked<1>(); // <1> is the number of array indices. See pybind11 docs.

    // 1D cell
    if (tdim == 1 && gdim == 1)
      {
        double p0 = mesh.coordinates()[v0];
        double p1 = mesh.coordinates()[v1];
        if (p0 > p1)
          std::swap(p0, p1);
        return p0 <= point[0] and point[0] <= p1;
      }

    // 2D cell
    else if (tdim == 2 && gdim == 2)
      {
        auto v0 = v(0);
        auto v1 = v(1);
        auto v2 = v(2);

        // mesh.coordinates() returns a 1D array of doubles: x0, y0, x1, y1,... in
        // this case.  So the value of x for vertex n is at offset n*2. Set pointers
        // to each vertex point:
        double* p0 = &mesh.coordinates()[v0*gdim];
        double* p1 = &mesh.coordinates()[v1*gdim];
        double* p2 = &mesh.coordinates()[v2*gdim];

        // The rest is code copied from CollisionPredicates.cpp
  
        const double ref = dnt::_orient2d(p0, p1, p2);

        if (ref > 0.0)
          {
            return (dnt::_orient2d(p1, p2, point) >= 0.0 and
                    dnt::_orient2d(p2, p0, point) >= 0.0 and
                    dnt::_orient2d(p0, p1, point) >= 0.0);
          }
        else if (ref < 0.0)
          {
            return (dnt::_orient2d(p1, p2, point) <= 0.0 and
                    dnt::_orient2d(p2, p0, point) <= 0.0 and
                    dnt::_orient2d(p0, p1, point) <= 0.0);
          }
        else
          {
            return ((dnt::_orient2d(p0, p1, point) == 0.0 and
                     dnt::_collides_segment_point_1d(p0[0], p1[0], point[0]) and
                     dnt::_collides_segment_point_1d(p0[1], p1[1], point[1])) or
                    (dnt::_orient2d(p1, p2, point) == 0.0 and
                     dnt::_collides_segment_point_1d(p1[0], p2[0], point[0]) and
                     dnt::_collides_segment_point_1d(p1[1], p2[1], point[1])) or
                    (dnt::_orient2d(p2, p0, point) == 0.0 and
                     dnt::_collides_segment_point_1d(p2[0], p0[0], point[0]) and
                     dnt::_collides_segment_point_1d(p2[1], p0[1], point[1])));
          }
      }

    // 3D cell
    else if (tdim == 3 && gdim == 3)
      {
        auto v0 = v(0);
        auto v1 = v(1);
        auto v2 = v(2);
        auto v3 = v(3);

        //  std::cout << "v0: " << v0 << " v1: " << v1 << " v2: " << v2 << " v3: " << v3 << std::endl;
        //  std::cout << "x: " << x << " y: " << y << " z: " << z << std::endl;
  
        // mesh.coordinates() returns a 1D array of doubles: x0, y0, z0, x1, y1, z1,... in this
        // case.  So the value of x for vertex n is at n*2.
        double* p0 = &mesh.coordinates()[v0*gdim];
        double* p1 = &mesh.coordinates()[v1*gdim];
        double* p2 = &mesh.coordinates()[v2*gdim];
        double* p3 = &mesh.coordinates()[v3*gdim];

        // std::cout << "p0: " << p0[0] << ", " << p0[1] << ", " << p0[2] << "\n"
        //           << "p1: " << p1[0] << ", " << p1[1] << ", " << p1[2] << "\n"
        //           << "p2: " << p2[0] << ", " << p2[1] << ", " << p2[2] << "\n"
        //           << "p3: " << p3[0] << ", " << p3[1] << ", " << p3[2]
        //           << std::endl;

        //  std::cout << "Point coords p: " << point[0] << " py: " << point[1] <<" pz: " << point[2] << std::endl;

        const double ref = dnt::_orient3d(p0, p1, p2, p3);

        if (ref > 0.0)
          {
            return (dnt::_orient3d(p0, p1, p2, point) >= 0.0 and
                    dnt::_orient3d(p0, p3, p1, point) >= 0.0 and
                    dnt::_orient3d(p0, p2, p3, point) >= 0.0 and
                    dnt::_orient3d(p1, p3, p2, point) >= 0.0);
          }
        else if (ref < 0.0)
          {
            return (dnt::_orient3d(p0, p1, p2, point) <= 0.0 and
                    dnt::_orient3d(p0, p3, p1, point) <= 0.0 and
                    dnt::_orient3d(p0, p2, p3, point) <= 0.0 and
                    dnt::_orient3d(p1, p3, p2, point) <= 0.0);
          }
        else
          {
            dolfin::dolfin_error("cell_contains_point (3D)",
                                 "compute tetrahedron point collision",
                                 "Not implemented for degenerate tetrahedron");
          }

        return false;
      }
  }    
  
  // 3D version using a DnT_pstruct
  // Test version, for a 1D mesh only

  bool cell_contains_point(dolfin::Mesh& mesh,
                           py::array_t<int> vertices,
                           Pstruct<Ptype::cartesian_xyz> point)
  {
    auto v = vertices.unchecked<1>(); // <1> is the number of array indices. See pybind11 docs.
  
    auto v0 = v(0);
    auto v1 = v(1);

    //  std::cout << "v0: " << v0 << " v1: " << v1 << std::endl;

    // 1D mesh version:
  
    double p0 = mesh.coordinates()[v0];
    double p1 = mesh.coordinates()[v1];

    // Copy the struct coordinates to a double array
    //  double point[] = {0.0, 0.0, 0.0};
    //  pstruct_to_point<Ptype::cartesian_xyz>(ps3D, point);

    if (p0 > p1)
      std::swap(p0, p1);
    return p0 <= point.x_ and point.x_ <= p1;
  }


  /// Versions that pass x, y, z directly as doubles

  // 1D based on passing x as a double
  bool cell_contains_point(dolfin::Mesh& mesh,
                           py::array_t<int> vertices,
                           double& x)
  {
    auto v = vertices.unchecked<1>(); // <1> is the number of array indices. See pybind11 docs.
  
    auto v0 = v(0);
    auto v1 = v(1);

    //  std::cout << "v0: " << v0 << " v1: " << v1 << std::endl;
    //  std::cout << "x: " << x << std::endl;
  
    double p0 = mesh.coordinates()[v0];
    double p1 = mesh.coordinates()[v1];

    if (p0 > p1)
      std::swap(p0, p1);
    return p0 <= x and x <= p1;
  }

  // 2D version based on passing x and y as doubles
  // Based on dolfin/dolfin/geometry/CollisionPredicates.cpp: collides_triangle_point_2d()

  bool cell_contains_point(dolfin::Mesh& mesh,
                           py::array_t<int> vertices,
                           double& x, double& y)
  {
    auto v = vertices.unchecked<1>(); // <1> is the number of array indices. See pybind11 docs.

    auto v0 = v(0);
    auto v1 = v(1);
    auto v2 = v(2);

    //  std::cout << "v0: " << v0 << " v1: " << v1 << " v2: " << v2 << std::endl;
    //  std::cout << "x: " << x << " y: " << y << std::endl;
  
    int meshDimension = 2;
    // mesh.coordinates() returns a 1D array of doubles: x0, y0, x1, y1,... in this
    // case.  So the value of x for vertex n is at n*2.
    double* p0 = &mesh.coordinates()[v0*meshDimension];
    double* p1 = &mesh.coordinates()[v1*meshDimension];
    double* p2 = &mesh.coordinates()[v2*meshDimension];

    // std::cout << "p0: " << p0[0] << ", " << p0[1] << "\n"
    //           << "p1: " << p1[0] << ", " << p1[1] << "\n"
    //           << "p2: " << p2[0] << ", " << p2[1]
    //           << std::endl;
  
    // Copy the struct coordinates to a double array
    double point[] = {x, y};
    //  std::cout << "p: " << point[0] << " py: " << point[1] << std::endl;

    // The rest is code copied from CollisionPredicates.cpp

    const double ref = dnt::_orient2d(p0, p1, p2);

    if (ref > 0.0)
      {
        return (dnt::_orient2d(p1, p2, point) >= 0.0 and
                dnt::_orient2d(p2, p0, point) >= 0.0 and
                dnt::_orient2d(p0, p1, point) >= 0.0);
      }
    else if (ref < 0.0)
      {
        return (dnt::_orient2d(p1, p2, point) <= 0.0 and
                dnt::_orient2d(p2, p0, point) <= 0.0 and
                dnt::_orient2d(p0, p1, point) <= 0.0);
      }
    else
      {
        return ((dnt::_orient2d(p0, p1, point) == 0.0 and
                 dnt::_collides_segment_point_1d(p0[0], p1[0], point[0]) and
                 dnt::_collides_segment_point_1d(p0[1], p1[1], point[1])) or
                (dnt::_orient2d(p1, p2, point) == 0.0 and
                 dnt::_collides_segment_point_1d(p1[0], p2[0], point[0]) and
                 dnt::_collides_segment_point_1d(p1[1], p2[1], point[1])) or
                (dnt::_orient2d(p2, p0, point) == 0.0 and
                 dnt::_collides_segment_point_1d(p2[0], p0[0], point[0]) and
                 dnt::_collides_segment_point_1d(p2[1], p0[1], point[1])));
      }

  }

  // 3D version based on passing x, y, z as doubles
  // Based on dolfin/dolfin/geometry/CollisionPredicates.cpp: collides_tetrahedron_point_3d()

  bool cell_contains_point(dolfin::Mesh& mesh,
                           py::array_t<int> vertices,
                           double& x, double& y, double& z)
  {
    auto v = vertices.unchecked<1>(); // <1> is the number of array indices. See pybind11 docs.

    auto v0 = v(0);
    auto v1 = v(1);
    auto v2 = v(2);
    auto v3 = v(3);

    //  std::cout << "v0: " << v0 << " v1: " << v1 << " v2: " << v2 << " v3: " << v3 << std::endl;
    //  std::cout << "x: " << x << " y: " << y << " z: " << z << std::endl;
  
    int meshDimension = 3;
    // mesh.coordinates() returns a 1D array of doubles: x0, y0, z0, x1, y1, z1,... in this
    // case.  So the value of x for vertex n is at n*2.
    double* p0 = &mesh.coordinates()[v0*meshDimension];
    double* p1 = &mesh.coordinates()[v1*meshDimension];
    double* p2 = &mesh.coordinates()[v2*meshDimension];
    double* p3 = &mesh.coordinates()[v3*meshDimension];

    // std::cout << "p0: " << p0[0] << ", " << p0[1] << ", " << p0[2] << "\n"
    //           << "p1: " << p1[0] << ", " << p1[1] << ", " << p1[2] << "\n"
    //           << "p2: " << p2[0] << ", " << p2[1] << ", " << p2[2] << "\n"
    //           << "p3: " << p3[0] << ", " << p3[1] << ", " << p3[2]
    //           << std::endl;
  
    // Copy the point coordinates to a double array
    double point[] = {x, y, z};
    //  std::cout << "Point coords p: " << point[0] << " py: " << point[1] <<" pz: " << point[2] << std::endl;

    const double ref = dnt::_orient3d(p0, p1, p2, p3);

    if (ref > 0.0)
      {
        return (dnt::_orient3d(p0, p1, p2, point) >= 0.0 and
                dnt::_orient3d(p0, p3, p1, point) >= 0.0 and
                dnt::_orient3d(p0, p2, p3, point) >= 0.0 and
                dnt::_orient3d(p1, p3, p2, point) >= 0.0);
      }
    else if (ref < 0.0)
      {
        return (dnt::_orient3d(p0, p1, p2, point) <= 0.0 and
                dnt::_orient3d(p0, p3, p1, point) <= 0.0 and
                dnt::_orient3d(p0, p2, p3, point) <= 0.0 and
                dnt::_orient3d(p1, p3, p2, point) <= 0.0);
      }
    else
      {
        dolfin::dolfin_error("cell_contains_point (3D)",
                             "compute tetrahedron point collision",
                             "Not implemented for degenerate tetrahedron");
      }

    return false;
  }

  //! Test whether a point lies within the cell given by its cell_index

  /*!
    Overload on the number of particle coordinates: 1D, 2D, 3D 
  */

  // 1D version Pstruct<Ptype::cartesian_x
  bool particle_stays_in_cell(Pstruct<Ptype::cartesian_x> p,
                              dolfin::Mesh& mesh,
                              py::dict cell_index_dict)
  {

    py::array_t vertices = cell_index_dict[p.cell_index_];
    auto v = vertices.unchecked<1>(); // <1> is the number of array indices. See pybind11 docs.
  
    auto v0 = v(0);
    auto v1 = v(1);

    //  std::cout << "v0: " << v0 << " v1: " << v1 << std::endl;
  
    double p0 = mesh.coordinates()[v0];
    double p1 = mesh.coordinates()[v1];
  
    if (p0 > p1)
      std::swap(p0, p1);
    return p0 <= p.x_ and p.x_ <= p1;
  }

  // 2D version Pstruct<Ptype::cartesian_xy
  // Based on dolfin/dolfin/geometry/CollisionPredicates.cpp: collides_triangle_point_2d()
  bool particle_stays_in_cell(Pstruct<Ptype::cartesian_xy> p,
                              dolfin::Mesh& mesh,
                              py::dict cell_index_dict)

  {
    py::array_t vertices = cell_index_dict[p.cell_index_];
    auto v = vertices.unchecked<1>(); // <1> is the number of array indices. See pybind11 docs.
  
    auto v0 = v(0);
    auto v1 = v(1);
    auto v2 = v(2);

    //  std::cout << "v0: " << v0 << " v1: " << v1 << std::endl;
  
    int meshDimension = 2;
    // mesh.coordinates() returns a 1D array of doubles: x0, y0, x1, y1,... in this
    // case.  So the value of x for vertex n is at n*2.
    double* p0 = &mesh.coordinates()[v0*meshDimension];
    double* p1 = &mesh.coordinates()[v1*meshDimension];
    double* p2 = &mesh.coordinates()[v2*meshDimension];

    // Copy the struct coordinates to a double array
    double point[] = {0.0, 0.0};
    pstruct_to_point<Ptype::cartesian_xy>(p, point);

    // The rest is code copied from CollisionPredicates.cpp
  
    const double ref = dnt::_orient2d(p0, p1, p2);

    if (ref > 0.0)
      {
        return (dnt::_orient2d(p1, p2, point) >= 0.0 and
                dnt::_orient2d(p2, p0, point) >= 0.0 and
                dnt::_orient2d(p0, p1, point) >= 0.0);
      }
    else if (ref < 0.0)
      {
        return (dnt::_orient2d(p1, p2, point) <= 0.0 and
                dnt::_orient2d(p2, p0, point) <= 0.0 and
                dnt::_orient2d(p0, p1, point) <= 0.0);
      }
    else
      {
        return ((dnt::_orient2d(p0, p1, point) == 0.0 and
                 dnt::_collides_segment_point_1d(p0[0], p1[0], point[0]) and
                 dnt::_collides_segment_point_1d(p0[1], p1[1], point[1])) or
                (dnt::_orient2d(p1, p2, point) == 0.0 and
                 dnt::_collides_segment_point_1d(p1[0], p2[0], point[0]) and
                 dnt::_collides_segment_point_1d(p1[1], p2[1], point[1])) or
                (dnt::_orient2d(p2, p0, point) == 0.0 and
                 dnt::_collides_segment_point_1d(p2[0], p0[0], point[0]) and
                 dnt::_collides_segment_point_1d(p2[1], p0[1], point[1])));
      }

  }

  // 3D version Pstruct<Ptype::cartesian_xyz
  // Based on dolfin/dolfin/geometry/CollisionPredicates.cpp: collides_tetrahedron_point_3d()
  bool particle_stays_in_cell(Pstruct<Ptype::cartesian_xyz> p,
                              dolfin::Mesh& mesh,
                              py::dict cell_index_dict)
  {
    py::array_t vertices = cell_index_dict[p.cell_index_];
    auto v = vertices.unchecked<1>(); // <1> is the number of array indices. See pybind11 docs.

    auto v0 = v(0);
    auto v1 = v(1);
    auto v2 = v(2);
    auto v3 = v(3);

    //  std::cout << "v0: " << v0 << " v1: " << v1 << " v2: " << v2 << " v3: " << v3 << std::endl;
    //  std::cout << "x: " << x << " y: " << y << " z: " << z << std::endl;
  
    int meshDimension = 3;
    // mesh.coordinates() returns a 1D array of doubles: x0, y0, z0, x1, y1, z1,... in this
    // case.  So the value of x for vertex n is at n*2.
    double* p0 = &mesh.coordinates()[v0*meshDimension];
    double* p1 = &mesh.coordinates()[v1*meshDimension];
    double* p2 = &mesh.coordinates()[v2*meshDimension];
    double* p3 = &mesh.coordinates()[v3*meshDimension];

    // std::cout << "p0: " << p0[0] << ", " << p0[1] << ", " << p0[2] << "\n"
    //           << "p1: " << p1[0] << ", " << p1[1] << ", " << p1[2] << "\n"
    //           << "p2: " << p2[0] << ", " << p2[1] << ", " << p2[2] << "\n"
    //           << "p3: " << p3[0] << ", " << p3[1] << ", " << p3[2]
    //           << std::endl;
  
    // Copy the struct coordinates to a double array
    double point[] = {0.0, 0.0, 0.0};
    pstruct_to_point<Ptype::cartesian_xy>(p, point);

    // The rest is code copied from CollisionPredicates.cpp
  
    const double ref = dnt::_orient3d(p0, p1, p2, p3);

    if (ref > 0.0)
      {
        return (dnt::_orient3d(p0, p1, p2, point) >= 0.0 and
                dnt::_orient3d(p0, p3, p1, point) >= 0.0 and
                dnt::_orient3d(p0, p2, p3, point) >= 0.0 and
                dnt::_orient3d(p1, p3, p2, point) >= 0.0);
      }
    else if (ref < 0.0)
      {
        return (dnt::_orient3d(p0, p1, p2, point) <= 0.0 and
                dnt::_orient3d(p0, p3, p1, point) <= 0.0 and
                dnt::_orient3d(p0, p2, p3, point) <= 0.0 and
                dnt::_orient3d(p1, p3, p2, point) <= 0.0);
      }
    else
      {
        dolfin::dolfin_error("particle_stays_in_cell (3D)",
                             "compute tetrahedron point collision",
                             "Not implemented for degenerate tetrahedron");
      }

    return false;
  }

  // We moved these to dolfin_cpp.cpp

  /*
    void dolfin_DnT_pstruct_instances()
    {
    dolfin::Mesh mesh;
    dolfin::Mesh& meshRef = mesh;
    py::array_t<int> vertices;
    Pstruct<Ptype::cartesian_x> ps1D;
    Pstruct<Ptype::cartesian_xy> ps2D;
    Pstruct<Ptype::cartesian_xyz> ps3D;

    // c++filt said that the following symbols were not in the .so file.  These
    // statements force the compiler to make these specialized versions of the
    // templated function in dolfin.o.

    // NEEDED?  
    cell_contains_point(meshRef, vertices, ps1D);
    cell_contains_point(meshRef, vertices, ps2D);
    cell_contains_point(meshRef, vertices, ps3D);
    }
  */

  //! Find the cell facet crossed in traveling along a displacement vector dr from position r0

  /*!
    Overload on the number of particle coordinates: 1D, 2D, 3D 
  */

  /// Versions using DnT_struct

  // 1D version using a DnT_pstruct
  //bool cell_contains_point(dolfin::Mesh& mesh,
  //                         py::array_t<int> vertices,
  //                         Pstruct<Ptype::cartesian_x> ps1D)

} // namespace dnt
