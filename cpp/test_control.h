// Copyright (C) 2019 L. D. Hughes

// Tests for pybind11 in C++

/*
Contents:

  function_with_DTcontrol_arg()

//remove:
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

*/

#ifndef TEST_CONTROL_H
#define TEST_CONTROL_H

namespace py = pybind11;

namespace dnt
{
  //! Test passing of various DnT Python classes to C++
  void function_with_DTcontrol_arg(py::object ctrl)
  {
    std::cout << "Hello from function_with_DTcontrol_arg" << std::endl;
    auto author = ctrl.attr("author").cast<std::string>();
    auto dt = ctrl.attr("dt").cast<double>();
    auto dt_str = std::to_string(dt);
    auto n_timesteps = ctrl.attr("n_timesteps").cast<int>();
    auto n_timesteps_str = std::to_string(n_timesteps);
    
    auto msg = "The author is " + author;
    std::cout << msg << std::endl;

    msg = "The timestep is " + dt_str;
    std::cout << msg << std::endl;

    msg = "The number of timesteps is " + n_timesteps_str;
    std::cout << msg << std::endl;
  }

} // namespace dnt

#endif
