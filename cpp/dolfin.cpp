#include "dolfin.h"
// #include "pstruct.h"

namespace py = pybind11;

/*
void init_ex2(py::module &m) {
m.def("sub", [](int a, int b) { return a - b; });
}
*/


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
template <typename PS, typename FS>
void interpolate_field_to_points(dolfin::Function& field,
                                 py::array_t<PS, 0> points,
                                 py::array_t<FS, 0> field_at_points)
{
  const auto points_info = points.request(); // request() returns metadata about the Python array (ptr, ndim, size, shape)
  const auto p = static_cast<PS*>(points_info.ptr); // Pointer to a struct in points
  
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

template <typename PS, typename FS>
void interpolate_field_to_points_v1(dolfin::Function& field,
                                 py::array_t<PS, 0> points,
                                 py::array_t<FS, 0> field_at_points)
{
  const auto points_info = points.request(); // request() returns metadata about the Python array (ptr, ndim, size, shape)
  const auto p = static_cast<PS*>(points_info.ptr); // Pointer to a struct in points
  
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
    pstruct_to_double<PS>(p[ip], point);
    
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

//! Test whether a point lies within the cell defined by the vertices.

bool cell_contains_point_1d(dolfin::Mesh& mesh,
                               py::array_t<int> vertices,
                               DnT_pstruct1D point)
{
  auto v = vertices.unchecked<1>(); // See pybind11 docs
  
  auto v0 = v(0);
  auto v1 = v(1);

//  std::cout << "v0: " << v0 << " v1: " << v1 << std::endl;
  
  double p0 = mesh.coordinates()[v0];
  double p1 = mesh.coordinates()[v1];
  
  if (p0 > p1)
    std::swap(p0, p1);
  return p0 <= point.x_ and point.x_ <= p1;
}

// Older versions:

bool cell_contains_point_1d_v1(dolfin::Mesh& mesh, unsigned int v0, unsigned int v1, double point)
{
//  std::size_t v0 = vertices[0];
//  std::size_t v1 = vertices[1];
  
// Use this for a Python list of ints called vertices[]:
  
//  auto v0 = vertices[0].cast<int>();
//  auto v1 = vertices[1].cast<int>();
  
//  auto v0 = vertices[0];
//  auto v1 = vertices[1];

//  cout << "v0: " << v0 << "v1: " << v1 << std::endl;
  
  double p0 = mesh.coordinates()[v0];
  double p1 = mesh.coordinates()[v1];
  
  if (p0 > p1)
    std::swap(p0, p1);
  return p0 <= point and point <= p1;
}

bool cell_contains_point_1d_v2(dolfin::Mesh& mesh,
                               py::list vertices,
                               DnT_pstruct1D point)

{
//  std::size_t v0 = vertices[0];
//  std::size_t v1 = vertices[1];
  
// Use this for a Python list of ints called vertices[]:
  auto v0 = vertices[0].cast<int>();
  auto v1 = vertices[1].cast<int>();
  
//  auto v0 = vertices[0];
//  auto v1 = vertices[1];

  std::cout << "v0: " << v0 << "v1: " << v1 << std::endl;
  
//  m.def("f_simple", [](DnT_pstruct1D p) { return p.x_ * 10; });


  
  double p0 = mesh.coordinates()[v0];
  double p1 = mesh.coordinates()[v1];
  
  if (p0 > p1)
    std::swap(p0, p1);
  return p0 <= point.x_ and point.x_ <= p1;
}

