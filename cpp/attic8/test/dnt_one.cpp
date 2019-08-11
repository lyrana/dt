//#include "dolfin.h"
//#include "pseg.h"

//namespace py = pybind11;


#include "pseg.h"


/*
void init_ex2(py::module &m) {
m.def("sub", [](int a, int b) { return a - b; });
}
*/

//! Sum the weights of a segment of particles in the cells they occupy.

/*! 

    This is the C++ version of the function add_weight_to_cell() called in
    Particle_Module.py. The function is templated on the type of the struct that
    corresponds to the Numpy pseg array. It loops over the particles in pseg and adds
    their weights to a cell-density array.

    \param pseg is a structured array containing particle data.

    \param dF is reference to a dolfin Function.

*/
template <typename PS>
void add_weights_to_cells(py::array_t<PS, 0> pseg, dolfin::Function& dF) {

  const auto pseg_info = pseg.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
  const auto p = static_cast<PS*>(pseg_info.ptr); // Pointer to a particle struct in pseg
  
  std::shared_ptr< const dolfin::FunctionSpace > dFS = dF.function_space();
  std::size_t dofs_per_cell = dFS->element()->space_dimension();

// This vector hold the values contributed by a particle to the density array. In the
// case of the cell density, a particle contributes to the one and only DoF in the
// cell, so we use the name 'weight' instead of 'weights':
  std::vector<double> weight(dofs_per_cell);
  
  std::shared_ptr< const dolfin::GenericDofMap > dDofmap = dFS->dofmap();

  // Print the DoF map. In this case it's cell indices to DoF indices.
  // dFS->print_dofmap();

  // Add the particle weights to the Vector
  for (auto i = 0; i < pseg_info.size; i++) {
    // std::cout << i << " x=" << p[i].x_ << std::endl;
    auto cellIndex = p[i].cell_index_;
    // std::cout << i << " cellIndex=" << cellIndex << std::endl;

    // In the case of a cell density, there's only one DoF in the cell, so we use the
    // name 'dofIndex' instead of 'dofIndices':
    auto dofIndex = dDofmap->cell_dofs(cellIndex);
    
    weight[0] = p[i].weight_;  // Only 1 DoF, so there's only 1 element in 'weight[]'.
    
    //dF.vector() is a std::shared_ptr<GenericVector>
    dF.vector()->add_local(weight.data(), dofs_per_cell, dofIndex.data());

// Q: would it be better to use the loop to fill arrays and then call this function
// only once, after the loop?

  }
  
}

// Define the C++ function interpolate_weights_to_dofs(): it's templated on the type
// of the struct that corresponds to the Numpy pseg array. It loops over each
// particles in pseg and adds its contribution to the DoFs in the cell contraining
// the particle.
template <typename PS>
void interpolate_weights_to_dofs(py::array_t<PS, 0> pseg, dolfin::Function& dF) { // The 0 means a continguous array with C ordering

  const auto pseg_info = pseg.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
  const auto p = static_cast<PS*>(pseg_info.ptr); // Pointer to a particle struct in pseg
  
  std::shared_ptr< const dolfin::FunctionSpace > dFS = dF.function_space();
  
  // Variables for cell information
  const dolfin::Mesh& mesh = *dFS->mesh();
  std::vector<double> dof_coordinates;
  ufc::cell ufc_cell;

  // Variables for evaluating basis
  const std::size_t rank = dFS->element()->value_rank();
  std::size_t size_basis = 1;
  for (std::size_t i = 0; i < rank; ++i)
    size_basis *= dFS->element()->value_dimension(i);
  std::size_t dofs_per_cell = dFS->element()->space_dimension();
  std::vector<double> basis(size_basis);

// This vector holds the values contributed by a particle to the density array.
  std::vector<double> weights(dofs_per_cell);
  
  // Variables for adding local information to vector
  double basis_sum;
  
  std::shared_ptr< const dolfin::GenericDofMap > dDofmap = dFS->dofmap();

  // Print the DoF map. In this case it maps cell indices to DoF indices.
  // dFS->print_dofmap();

  // Add the particle weights to the Vector
  double point[] = {0.0, 0.0, 0.0};
  for (auto ip = 0; ip < pseg_info.size; ip++) {
    // std::cout << ip << " x=" << p[ip].x_ << std::endl;

    // Copy spatial coordinates of p[ip] to a double[] for passing to evaluate_basis()
    pstruct_to_double<PS>(p[ip], point);
    
    auto cellIndex = p[ip].cell_index_;
    // std::cout << ip << " cellIndex=" << cellIndex << std::endl;

    // Get the indices of the DoFs in this cell.
    auto dofIndices = dDofmap->cell_dofs(cellIndex);
    
    // Get the cell object containing this particle
    dolfin::Cell cell(mesh, static_cast<std::size_t>(cellIndex));
    // Get the coordinates of the DoFs in this cell
    cell.get_coordinate_dofs(dof_coordinates);

    // Evaluate all basis functions at the point[]
    cell.get_cell_data(ufc_cell);
    for (std::size_t i = 0; i < dofs_per_cell; ++i)
    {
      dFS->element()->evaluate_basis(i, basis.data(),
                                     point,  // const double* x
                                     dof_coordinates.data(),
                                     ufc_cell.orientation);
      basis_sum = 0.0;
      for (const auto& v : basis)
        basis_sum += v;
      weights[i] = p[ip].weight_*basis_sum;
    }

    //dF.vector() is a std::shared_ptr<GenericVector>
    dF.vector()->add_local(weights.data(), dofs_per_cell, dofIndices.data());
    
  }
  
}











// The PYBIND11_MODULE() macro creates a function that will be called when an import
// statement is issued from within Python. The module name (dnt_one) is given as the
// first macro argument (it should not be in quotes). The second argument (m) defines
// a variable of type py::module which is the main interface for creating
// bindings. The method module::def() generates binding code that exposes the C++
// functions to Python.

// // Create a variable 'm' of type py::module
PYBIND11_MODULE(dnt_one, m) {

// typeinfo may be registered before the dtype descriptor for scalar casts to work...


///
// Allow Python types to be used in C++
///  
  
//  py::class_() creates Python bindings for a C++ class or struct-style data
// structure.  The following statements connect the Python names DnT_pstructXD to C++
// structs with the same names.
  py::class_<DnT_pstruct1D>(m, "DnT_pstruct1D");
//  py::class_<DnT_pstruct2D>(m, "DnT_pstruct2D");
//  py::class_<DnT_pstruct3D>(m, "DnT_pstruct3D");

// The following statements allow the Numpy structured types DnT_pstructXD to be used in C++.
// In particular they can be used as template arguments to py::array_t.
  
// Register DnT_pstruct1D as a Numpy dtype descriptor. DnT_pstruct1D is of type py::dtype (a class).
//  PYBIND11_NUMPY_DTYPE(DnT_pstruct1D, x, x0, ux, weight, bitflags, cell_index, unique_ID, crossings);
// The "_EX" variation allows the Python names of the variables in the structure to be different from the variable names in the C++ struct.
  PYBIND11_NUMPY_DTYPE_EX(DnT_pstruct1D, x_, "x", x0_, "x0", ux_, "ux", weight_, "weight", bitflags_, "bitflags", cell_index_, "cell_index", unique_ID_, "unique_ID", crossings_, "crossings");

// Register DnT_pstruct2D and DnT_pstruct3D
//  PYBIND11_NUMPY_DTYPE_EX(DnT_pstruct2D, x_, "x", y_, "y", x0_, "x0", y0_, "y0", ux_, "ux", uy_, "uy", weight_, "weight", bitflags_, "bitflags", cell_index_, "cell_index", unique_ID_, "unique_ID", crossings_, "crossings");

//  PYBIND11_NUMPY_DTYPE_EX(DnT_pstruct3D, x_, "x", y_, "y", z_, "z", x0_, "x0", y0_, "y0", z0_, "z0", ux_, "ux", uy_, "uy", uz_, "uz", weight_, "weight", bitflags_, "bitflags", cell_index_, "cell_index", unique_ID_, "unique_ID", crossings_, "crossings");

  
///
// Bindings for particle-to-mesh operations
///  
  
// Connect the Python symbol add_weights_to_cells() to the C++ function declared as
//       void add_weights_to_cells(py::array_t<S, 0> pseg, dolfin::Function& dF) {}
//  m.def("add_weights_to_cells1D", (void * (py::array_t<DnT_pstruct1D, 0>, dolfin::Function&)) &add_weights_to_cells<DnT_pstruct1D>);
  
//  m.def("add_weights_to_cells1D", (void (*) (pybind11::array_t<DnT_pstruct1D, 0>, dolfin::Function&)) &add_weights_to_cells<DnT_pstruct1D>);
  m.def("add_weights_to_cells1D",  &add_weights_to_cells<DnT_pstruct1D>);

  
//  m.def("add_weights_to_cells2D", &add_weights_to_cells<DnT_pstruct2D>);
//  m.def("add_weights_to_cells3D", &add_weights_to_cells<DnT_pstruct3D>);

// Connect the Python symbol interpolate_weights_to_dofs to the C++ function declared as
//       void interpolate_weights_to_dofs(py::array_t<S, 0> pseg, dolfin::Function& dF) {}
//  m.def("interpolate_weights_to_dofs1D", &interpolate_weights_to_dofs<DnT_pstruct1D>);
//  m.def("interpolate_weights_to_dofs2D", &interpolate_weights_to_dofs<DnT_pstruct2D>);
//  m.def("interpolate_weights_to_dofs3D", &interpolate_weights_to_dofs<DnT_pstruct3D>);

  
///
// Bindings for mesh-to-particle functions
///


///
// Bindings for printing functions
///
  
// Connect the Python symbol print_pseg() to the C++ function declared as:
//                   template <typename S>
//                   py::list print_pstructarray(py::array_t<S, 0> arr) {}
//  m.def("print_pseg1D", &print_pstructarray<DnT_pstruct1D>);
//  m.def("print_pseg2D", &print_pstructarray<DnT_pstruct2D>);
//  m.def("print_pseg3D", &print_pstructarray<DnT_pstruct3D>);

}

