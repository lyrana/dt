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

  // Print the DoF map. In this case it's a map of cell indices to DoF indices.
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
// particle in pseg and adds its contribution to the DoFs in the cell contraining
// the particle.
template <typename PS>
void interpolate_weights_to_dofs(py::array_t<PS, 0> pseg, dolfin::Function& dF) { // The 0 means a continguous array with C ordering

  // NEED to pass in npSeg to limit the range of the particle loop !!!!!!
  // OR, pass pseg[0:npSeg] instead?

  const auto pseg_info = pseg.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
  const auto p = static_cast<PS*>(pseg_info.ptr); // Pointer to a particle struct in pseg

  std::shared_ptr< const dolfin::FunctionSpace > dFS = dF.function_space();

  // Variables for cell information
  const dolfin::Mesh& mesh = *dFS->mesh();
  std::vector<double> dof_coordinates;
  ufc::cell ufc_cell;

  // Variables for evaluating the cell basis functions
  const std::size_t rank = dFS->element()->value_rank();
  std::size_t size_basis = 1;
  for (std::size_t i = 0; i < rank; ++i)
    size_basis *= dFS->element()->value_dimension(i);
  std::size_t dofs_per_cell = dFS->element()->space_dimension();
  std::vector<double> basis(size_basis);

// This vector holds the values contributed by a particle to the density array.
  std::vector<double> weights(dofs_per_cell);
  
  std::shared_ptr< const dolfin::GenericDofMap > dDofmap = dFS->dofmap();

  // Print the DoF map. In this case it maps cell indices to DoF indices.
  // dFS->print_dofmap();

  // Add the particle weights to the Vector
  double point[] = {0.0, 0.0, 0.0};
  for (auto ip = 0; ip < pseg_info.size; ip++) {
    // std::cout << ip << " x=" << p[ip].x_ << std::endl;

    // Copy spatial coordinates of p[ip] to a double[] for passing to evaluate_basis()

    //TODO: Could p[] be used directly, avoiding the copy?
    pstruct_to_point<PS>(p[ip], point);
    
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
      double basis_sum = 0.0;
      for (const auto& v : basis)
        basis_sum += v;
      weights[i] = p[ip].weight_*basis_sum;
    }

    //dF.vector() is a std::shared_ptr<GenericVector>
    dF.vector()->add_local(weights.data(), dofs_per_cell, dofIndices.data());
    
  }
  
}

//! Push a segment of neutral particles

/*!
    Overload on the number of particle coordinates: 1D, 2D, 3D 
 */

// or

/*

template <typename PS>
  void move_neutral_particle_segment(py::array_t<PS, 0> psegIn,py::array_t<PS, 0> psegOut, dolfin::Function& dF) { // The 0 means a continguous array with C ordering

  const auto pseg_info = pseg.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
  const auto p = static_cast<PS*>(pseg_info.ptr); // Pointer to a particle struct in pseg
  
  std::shared_ptr< const dolfin::FunctionSpace > dFS = dF.function_space();
  
  // Variables for cell information
  const dolfin::Mesh& mesh = *dFS->mesh();
  std::vector<double> dof_coordinates;
  ufc::cell ufc_cell;

  // Variables for evaluating the cell basis functions
  const std::size_t rank = dFS->element()->value_rank();
  std::size_t size_basis = 1;
  for (std::size_t i = 0; i < rank; ++i)
    size_basis *= dFS->element()->value_dimension(i);
  std::size_t dofs_per_cell = dFS->element()->space_dimension();
  std::vector<double> basis(size_basis);

// This vector holds the values contributed by a particle to the density array.
  std::vector<double> weights(dofs_per_cell);
  
  std::shared_ptr< const dolfin::GenericDofMap > dDofmap = dFS->dofmap();

  // Print the DoF map. In this case it maps cell indices to DoF indices.
  // dFS->print_dofmap();

  // Add the particle weights to the Vector

  
  //// The particle loop

  
  double point[] = {0.0, 0.0, 0.0};
  for (auto ip = 0; ip < pseg_info.size; ip++) {
    // std::cout << ip << " x=" << p[ip].x_ << std::endl;

    // Copy spatial coordinates of p[ip] to a double[] for passing to evaluate_basis()

    //TODO: Could p[] be used directly, avoiding the copy?
    pstruct_to_point<PS>(p[ip], point);
    
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
      double basis_sum = 0.0;
      for (const auto& v : basis)
        basis_sum += v;
      weights[i] = p[ip].weight_*basis_sum;
    }

    //dF.vector() is a std::shared_ptr<GenericVector>
    dF.vector()->add_local(weights.data(), dofs_per_cell, dofIndices.data());
    
  }
  
}

*/

void pseg_DnT_pstruct_instances()
{
  py::array_t<DnT_pstruct1D, 0> pseg1D;
  py::array_t<DnT_pstruct2D, 0> pseg2D;
  py::array_t<DnT_pstruct3D, 0> pseg3D;
  dolfin::Function dolfinFunction;
  dolfin::Function& dFref = dolfinFunction;

  // c++filt said that the following symbols were not in the .so file.  These
  // statements force the compiler to make these specialized versions of the
  // templated function in pseg.o.
  add_weights_to_cells<DnT_pstruct1D>(pseg1D, dFref);
  add_weights_to_cells<DnT_pstruct2D>(pseg2D, dFref);
  add_weights_to_cells<DnT_pstruct3D>(pseg3D, dFref);
  
  interpolate_weights_to_dofs<DnT_pstruct1D>(pseg1D, dFref);
  interpolate_weights_to_dofs<DnT_pstruct2D>(pseg2D, dFref);
  interpolate_weights_to_dofs<DnT_pstruct3D>(pseg3D, dFref);
}
