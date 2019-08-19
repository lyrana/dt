// Copyright (C) 2019 L. D. Hughes

// Utility functions for DnT in C++

/*
Contents:

  dnt_error()

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

#ifndef DNT_H
#define DNT_H

namespace dnt

{

  // Based on dolfin/log/log/Logger.cpp
  //void Logger::error(std::string msg) const  
  //! Handle a fatal error
  void dnt_error(std::string msg) const
  {
    std::string s = std::string("DnT error message: ") + msg;
    throw std::runtime_error(s);
  }

} // namespace dnt

#endif
