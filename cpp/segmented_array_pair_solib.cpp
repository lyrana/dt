/*! \file segmented_array_pair_solib.cpp

  \brief This file creates a shared library containing the Python bindings for
  SegmentedArrayPair storage that's created in C++.

  \namespace dnt

  This file contains the Python-to-C++ bindings allows access to the C++ version of
  SegmentedArrayPair storage.

*/
#include "SegmentedArrayPair.h"

namespace py = pybind11;

// The PYBIND11_MODULE() macro creates a function that will be called when an import
// statement is issued from within Python. The module name (dnt_cpp) is given as the
// first macro argument (it should not be in quotes). The second argument (m) defines
// a variable of type py::module which is the main interface for creating
// bindings. The method module::def() generates binding code that exposes the C++
// functions to Python.


namespace dnt {

  // The anonymous namespace limits the scope of the functions in it to this file.
  namespace {

    //! Make SegmentedArray storage for particular particle types.
    /*!

      makeSegmentedArrayPair() is a 'helper' function. It causes the compiler to
      make the specialized SegmentedArrayPair storage class and member functions for the
      type of particle struct specified by a Ptype template parameter. These are then
      callable from Python.

      \param PT specifies the Ptype of the particle struct that needs to be stored.
      \param m is a py::module object created by PYBIND11_MODULE
      \param PStype is a string used to identify the particle type being stored.
      \return void
      \sa Ptype, SegmentedArrayPair

     */
    template <Ptype PT>
    void makeSegmentedArrayPair(py::module &m, std::string const & PStype)
    {
      using SAP = SegmentedArrayPair<PT>;
      // Make a class name with the particle structure type appended to
      // "SegmentedArrayPair"
      std::string pyclass_name = std::string("SegmentedArrayPair_") + PStype;
      
      // Create the Python binding for 
      py::class_<SAP>(m, pyclass_name.c_str())

        // The ctor
        
        // The C++ ctor args types are template parameters in py::init<>()
        // Here, the parameter is the segment length.
        .def(py::init<py::ssize_t>())

        // The following creates the bindings to SegmentedArrayPair member
        // functions for particle type PT. The source code for these is in
        // SegmentedArrayPair.h.
        .def("get_as_list", &SAP::get_as_list)
        .def("get_as_tuple", &SAP::get_as_tuple)
        .def("get_capacity", &SAP::get_capacity)
        .def("get_full_index", &SAP::get_full_index)
        .def("get_full_indices", &SAP::get_full_indices)
        .def("get_item", &SAP::get_item)
        // Set a default value of a keyword argument:
        .def("get_next_out_segment", &SAP::get_next_out_segment, py::arg("return_data_ptr") = false)
        .def("get_next_segment", &SAP::get_next_segment, py::arg(), py::arg("return_data_ptr") = false)
        .def("get_current_out_segment", &SAP::get_current_out_segment)
        .def("get_number_of_items", &SAP::get_number_of_items)
        .def("get_number_of_mbytes", &SAP::get_number_of_mbytes)
        .def("get_number_of_segments", &SAP::get_number_of_segments)
        .def("get_segment_and_offset", &SAP::get_segment_and_offset)
        .def("init_inout_loop", &SAP::init_inout_loop, py::arg("return_data_ptrs") = false)
        .def("init_out_loop", &SAP::init_out_loop, py::arg("return_data_ptr") = false)
        // Overloaded versions of push_back().
        // Put a signature in parens in front of &SAP::push_back to disambiguate the versions.
        // Note: (py::tuple (SAP::*) (arg)) &SAP::push_back(): This means that &SAP::push_back
        // is a pointer-to-a-function that takes 'arg' and returns a py::tuple.
        .def("push_back", (py::tuple (SAP::*) (py::array_t<Pstruct<PT>, 0>)) &SAP::push_back)
        .def("push_back", (py::tuple (SAP::*) (py::tuple)) &SAP::push_back)
        .def("push_back", (py::tuple (SAP::*) (py::list)) &SAP::push_back)
        .def("set_number_of_items", &SAP::set_number_of_items);
             
      // cf. ~/workspace/dolfin/python/src/geometry.cpp for __getitem__
        /*
        .def("__getitem__", [](SAP& self, std::size_t full_index)
             {
               //               if (index > 2)
               //                 throw py::index_error("Out of range");

               // (seg, offset) = divmod(i, self.SEGMENTLENGTH)
               auto segmentLength = self.get_segment_length();
               auto seg = (py::ssize_t) full_index / segmentLength;
               auto offset = full_index % segmentLength;
      
               auto outSA = self.outSegmentedArray;

               py::buffer_info pseg_info = self.segListPair[outSA][seg].request(); // request() returns metadata about the array (ptr, ndim, size, shape)
               const auto pseg = static_cast<Pstruct<PT>*>(pseg_info.ptr); // Pointer to a  structured Numpy array

               //        vec[firstAvailableOffset[outSA]] = item_input;
               //               auto ip = firstAvailableOffset[outSA];
               
               return self.segListPair[outSA][seg][offset];
               
               //               return self[index];
             });
        */
        
    } // void makeSegmentedArrayPair(py::module &m, std::string const & PStype)

  } // namespace none
  
  // Interface to the C++ class SegmentedArrayPair

  // Create a variable 'm' of type py::module
  PYBIND11_MODULE(segmented_array_pair_solib, m)
  {
    
    // The following registrations are needed for the Pstruct<> structured types to work with py::array_t. See pybind11: 12.2.3 Structured types.
    // The "_EX" variation allows the Python names of the variables in the structure to be different from the variable names in the C++ struct.
    PYBIND11_NUMPY_DTYPE_EX(Pstruct<Ptype::cartesian_x>, x_, "x", x0_, "x0", ux_, "ux", weight_, "weight", bitflags_, "bitflags", cell_index_, "cell_index", unique_ID_, "unique_ID", crossings_, "crossings");
    PYBIND11_NUMPY_DTYPE_EX(Pstruct<Ptype::cartesian_xy>, x_, "x", y_, "y", x0_, "x0", y0_, "y0", ux_, "ux", uy_, "uy", weight_, "weight", bitflags_, "bitflags", cell_index_, "cell_index", unique_ID_, "unique_ID", crossings_, "crossings");
    PYBIND11_NUMPY_DTYPE_EX(Pstruct<Ptype::cartesian_xyz>, x_, "x", y_, "y", z_, "z", x0_, "x0", y0_, "y0", z0_, "z0", ux_, "ux", uy_, "uy", uz_, "uz", weight_, "weight", bitflags_, "bitflags", cell_index_, "cell_index", unique_ID_, "unique_ID", crossings_, "crossings");

    // Create the specific classes "SegmentedArrayPair_cartesian_x/_xy/_xyz"
    makeSegmentedArrayPair<Ptype::cartesian_x>(m, "cartesian_x");
    makeSegmentedArrayPair<Ptype::cartesian_xy>(m, "cartesian_xy");
    makeSegmentedArrayPair<Ptype::cartesian_xyz>(m, "cartesian_xyz");
      
  } // ENDDEF: PYBIND11_MODULE(segmented_array_pair_solib, m)

  
} // namespace dnt
