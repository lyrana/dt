/*! \file particle.h

  \brief This file has the source code for a C++ implementation of particle advance.

  \namespace dnt

  template <Ptype PT>
    void move_neutral_particle_species(py::object species_name, py::object ctrl, py::object particle_P)


  \sa SegmentedArrayPair
*/

#ifndef PARTICLE_H
#define PARTICLE_H

#include <dolfin/common/Array.h>
#include <dolfin/common/Hierarchical.h>
#include <dolfin/function/assign.h>
#include <dolfin/function/Constant.h>
#include <dolfin/function/Expression.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionAssigner.h>
#include <dolfin/function/FunctionAXPY.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/MultiMeshFunction.h>
#include <dolfin/function/MultiMeshFunctionSpace.h>
#include <dolfin/function/MultiMeshSubSpace.h>
#include <dolfin/function/LagrangeInterpolator.h>
#include <dolfin/function/SpecialFunctions.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/fem/MultiMeshDofMap.h>
#include <dolfin/geometry/Point.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/mesh/Mesh.h>

#include <iostream>
#include <tuple>
//#include "pstruct.h"
#include "SegmentedArrayPair.h"

#include "dolfin.h"

//using namespace std; 

//template<typename PS>
//using SegList = std::vector<py::array_t<PS,0>>;

namespace py = pybind11;

namespace dnt
{
  //! Advance a neutral particle species by one time-increment on a mesh.
  /*!          

    Compute change in position in time dt. Use an explicit method to calculate the
    final position.  If a particle leaves its initial cell, the cell that the
    particle moves to is calculated by finding what facet the particle crosses, and
    looking up what the neighbor cell is. This is repeated until the cell
    containing the final position is found.

    \param ctrl is a DTcontrol_C object.
    \param psa is a SegmentedArrayPair containing the particle data for one neutral species.
    \param mesh is a dolfin Mesh that the particles move through.
    \cell_vertex_dict is a dictionary that returns the vertex ids, given a cell id.

    :cvar double dtRemaining: The time a particle has left to move in the new cell
                              the particle has entered.
    :cvar int pDim: Number of spatial coordinates in the particle location.
    :cvar double tStart: The time at which a particle starts its move in the
                                 current cell.

  */
  template <Ptype PT>
  //  void move_neutral_particle_species(py::str species_name, py::object ctrl, py::object particle_P)
  void move_neutral_particle_species(py::object species_name, py::object ctrl, py::object particle_P)
  {

    std::cout << "Hello from move_neutral_particle_species@0" << std::endl;

    // Get parameters needed from ctrl arg
    auto dt = ctrl.attr("dt").cast<double>(); // This does a COPY.
    auto step = ctrl.attr("timeloop_count").cast<int>();
    auto t = ctrl.attr("time").cast<double>();
    auto MAX_FACET_CROSS_COUNT = ctrl.attr("MAX_FACET_CROSS_COUNT").cast<int>();
    
    std::cout << "MAX_FACET_CROSS_COUNT " << MAX_FACET_CROSS_COUNT << std::endl;
    
    // Get attributes needed from Particle_C arg
    auto pmesh_M = particle_P.attr("pmesh_M");
    auto mesh = pmesh_M.attr("mesh").cast<dolfin::Mesh>();
    auto NO_FACET = pmesh_M.attr("NO_FACET").cast<int>();

    std::cout << "NO_FACET " << NO_FACET << std::endl;
    
    // py::array_t<int>* vertices2 = pmesh_M.attr("cell_vertex_dict")["2"]; // This doesn't compile.
    // auto vertices1 = pmesh_M.attr("cell_vertex_dict")["1"].cast<py::array_t<int>>(); // This compiles but get Keyerror on execute.
    //auto vertices2 = pmesh_M.attr("cell_vertex_dict")[2].cast<py::array_t<int>>(); // This doesn't compile.
    //auto v2 = vertices2.unchecked<1>();
    //std::cout << "v2 " << v2(0) << ", " << v2(1) << std::endl;
    
    //    py::print("vertices2: ", py::str(vertices2[2])); // Need to use the buffer protocol to get at the data.
    
    //v1: auto cell_vertex_dict = pmesh_M.attr("cell_vertex_dict").cast<std::map<unsigned int, py::array_t<int>>>(); // Does this do a copy op?

    // This isn't a correct way to get at the py::dict!
    //    auto cellVertexDictCap = pmesh_M.attr("cell_vertex_dict").cast<py::capsule>();    
    //    Eigen::Map<const Eigen::Matrix<unsigned int, Eigen::Dynamic, 1>>* cell_vertex_dict = cellVertexDictCap;
      
    //std::cout << "vertices " << cell_vertex_dict[2];
    //v1: py::print("cell_vertex_dict[2] ", py::str(cell_vertex_dict[2]));

    //    std::cout << "cell_vertex_dict " << cell_vertex_dict << std::endl;
    
    auto pDim = particle_P.attr("particle_dimension").cast<int>();

    std::cout << "pDim " << pDim << std::endl;

    //1 auto pseg_arr = particle_P.attr("pseg_arr").cast<py::dict>();

    std::cout << "Hello from move_neutral_particle_species@1" << std::endl;
    
    // Get the data for the species to be advanced
    auto pseg_arr_map = particle_P.attr("pseg_arr").cast<std::map<std::string, SegmentedArrayPair<PT> *>>();
    auto psa = pseg_arr_map[std::string(py::str(species_name))];

    // Access the SAP without constructing a std::map. Note that psa is a pointer.
    //1 auto psa = pseg_arr[species_name].cast<SegmentedArrayPair<PT> *>();

    const py::ssize_t segmentLength = psa->get_segment_length();
    //    auto mesh_coordinates = mesh.coordinates();
    
    // Start a loop over the psa segments. segTuple contains (npSeg, psegIn, psegOut)
    py::tuple segTuple = psa->init_inout_loop(true); // 'true' returns the data pointers for
                                                     // the first segment pair, instead of the
                                                     // Numpy arrays.
    auto npSeg = segTuple[0].cast<py::ssize_t>();

    // Get the particle array pointers from the returned py::capsule's
    // psegIn:
    auto psegInCap = segTuple[1].cast<py::capsule>();
    // std::cout << "The name of the capsule is " << psegInCap.name() << std::endl;
    Pstruct<PT>* psegIn = psegInCap; // Declares the type of psegIn
    // std::cout << "psegIn[0].x_ = " << psegIn[0].x_ << std::endl;
    // psegOut:
    auto psegOutCap = segTuple[2].cast<py::capsule>();
    // std::cout << "The name of the capsule is " << psegOutCap.name() << std::endl;
    Pstruct<PT>* psegOut = psegOutCap; // Declares the type of psegOut

    // An alternative way to do this is to use the returned Numpy structured array
    // objects: call init_inout_loop() without the "true" argument above. Then use the
    // Python buffer protocol to get at the data.
    /*
    // psegIn:
    auto psegInNumpyArray = segTuple[1].cast<py::array_t<Pstruct<PT>,0>>();
    const auto psegInNumpyArray_info = psegInNumpyArray.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
    auto psegIn = static_cast<Pstruct<PT>*>(psegInNumpyArray_info.ptr); // Pointer to a particle struct in psegIn
    // psegOut:
    auto psegOutNumpyArray = segTuple[2].cast<py::array_t<Pstruct<PT>,0>>();
    const auto psegOutNumpyArray_info = psegOutNumpyArray.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
    auto psegOut = static_cast<Pstruct<PT>*>(psegOutNumpyArray_info.ptr); // Pointer to a particle struct in psegOut
    */
    
    py::ssize_t particleCount = 0; // Counts the number of still-active particles for this
                                   // species.
    
    // ipOut counts particles being written to the current "out" segment.
    py::ssize_t ipOut = 0;
    bool indexChange = false; // Flag to indicate if the particle SA indices have
                              // changed. This affects, e.g., trajectories, which use SA
                              // indices to identify particles chosen for trajectory
                              // plots.

    std::cout << "Hello from move_neutral_particle_species@2" << std::endl;
    
    while (psegIn != nullptr) // Keep looping until we run out of "in" segments.
      {
        for (auto ipIn = 0; ipIn < npSeg; ipIn++)
          {
            // Skip deleted particles
            //            if (psegIn[ipIn]['bitflags'] & self.DELETE_FLAG != 0)
            if ((psegIn[ipIn].bitflags_ & Pstruct<PT>::DELETE_FLAG) != 0b0)
              {
                indexChange = true; // Particle SA indices are stale past this point
                                    // due to deletions.
                continue; // skip to the next particle
              }

            // If this "out" segment is full, get the next "out" segment, since
            // there are still some "in" particles to advance.  If there are no more
            // "out" segments, allocate a new one.
            if (ipOut == segmentLength)
              {
                // Using pointers to the plain arrays:
                py::tuple segTuple = psa->get_next_out_segment(true);                
                // auto psegOutCap = segTuple[0].cast<py::capsule>();
                // Pstruct<PT>* psegOut = psegInCap;
                psegOut = segTuple[0].cast<py::capsule>(); // This compresses the above 2 lines into 1
                
                // Using Numpy structured arrays instead:
                /*
                py::tuple segTuple = psa->get_next_out_segment();
                auto psegOutNumpyArray = segTuple[0].cast<py::array_t<Pstruct<PT>,0>>();
                const auto psegOutNumpyArray_info = psegOutNumpyArray.request();
                const auto psegOut = static_cast<Pstruct<PT>*>(psegOutNumpyArray_info.ptr);
                */
                
                ipOut = 0; // Reset the counter for the new segment
              }
            
            psegOut[ipOut] = psegIn[ipIn]; // Copy all of this particle's data from the
                                           // input slot to the output slot

            // Move the neutral particle
            // 1. Save the initial position
            
            //            std::cout << "" << std::endl;
            std::cout << "ipOut " << ipOut << std::endl;
            
            psegOut[ipOut].x0_ = psegOut[ipOut].x_;
            psegOut[ipOut].y0_ = psegOut[ipOut].y_;
            psegOut[ipOut].z0_ = psegOut[ipOut].z_;
            // 2. Move the particle. NON-RELATIVISTIC: v and u are the same
            psegOut[ipOut].x_ += psegOut[ipOut].ux_*dt;
            psegOut[ipOut].y_ += psegOut[ipOut].uy_*dt;
            psegOut[ipOut].z_ += psegOut[ipOut].uz_*dt;

            auto pCellIndex = psegOut[ipOut].cell_index_;

            std::cout << "pCellIndex " << pCellIndex << std::endl;            

            // Loop until the particle is in the current cell
            auto mLastFacet = NO_FACET;
            size_t facetCrossCount = 0;
            auto tStart = t - dt;
            auto dtRemaining = dt;

            /*            
            std::string sIndex = std::to_string(pCellIndex);
            char const* pchar = sIndex.c_str();
            //            auto vertices = cell_vertex_dict[pchar].cast<Eigen::Matrix>();

            std::cout << "pchar = " << pchar << std::endl;
            */

            /*            
            for (auto item : cell_vertex_dict)
              {
                std::cout << "key=" << std::string(py::str(item.first)) << ", "
                          << "value=" << std::string(py::str(item.second)) << std::endl;
              }
            */

            //v1: auto vertices = cell_vertex_dict[pCellIndex]; // this returns a py::array_t<int>
            // Newer version that calls directly to C++
            // Create a view to this cell in the dolfin::Mesh object
            dolfin::Cell pcell(mesh, pCellIndex);
            auto vertices = pcell.entities(0);

            // std::cout << "vertices " << vertices[0] << ", " << vertices[1] << std::endl;
            
            double point[3] = {psegOut[ipOut].x_, psegOut[ipOut].y_, psegOut[ipOut].z_};
            
            std::cout << "Calling cell_contains_point" << std::endl;
            
            while (! cell_contains_point(mesh, vertices, point))
              {
                std::cout << "Hello from call to cell_contains_point" << std::endl;
                // Copy code here from particle-v1.cpp
                
              }
            
            /*
              Check if we've reached the end of this segment.  If so, we need to start
              writing on a new segment.  If there are no more segments, allocate a new one.
            */
            if (ipOut == segmentLength)
              {
                particleCount += segmentLength;
              }

          } // loop over particles in the "in" segment
        
        // Done with this "in" segment. Get the next one, if it exists.
        py::tuple segTuple = psa->get_next_segment("in", true);
        auto npSeg = segTuple[0].cast<py::ssize_t>();
        if (npSeg >  0)
          {
            psegIn = segTuple[1].cast<py::capsule>(); // This is compressed: there's a cast performed to make psegIn, which has it's type from above.
          }
        else
          {
            // py::print("*segTuple: ", *segTuple);
            psegIn = segTuple[1].cast<nullptr_t>(); // This will fail if segTuple[1] is not None
            assert(psegIn == nullptr);
          }

        // Using Numpy structured arrays:
        /*
        py::tuple segTuple = psa->get_next_segment("in");
        npSeg = segTuple[0].cast<py::ssize_t>();
        auto psegInNumpyArray = segTuple[1].cast<py::array_t<Pstruct<PT>,0>>();
        const auto psegInNumpyArray_info = psegInNumpyArray.request();
        psegIn = static_cast<Pstruct<PT>*>(psegInNumpyArray_info.ptr);
        */
      } // while there are more "in" segments, keep looping

    // Set values that will be used to keep track of the particle arrays
    if (ipOut != segmentLength) // This catches the case where we exit the loop when ipOut = segmentLength and there are no more "in" segments.  Otherwise, we would add segmentLength to particleCount twice.
      {
        particleCount += ipOut;
      }
    psa->set_number_of_items("out", particleCount);
  }
  // ENDDEF: void move_neutral_particle_species(py::object species_name, py::object ctrl, py::object particle_P)
  
  //! Advance a neutral particle species by one time increment on a mesh.
  /*!          

    Compute change in position in time dt. Use an explicit method to calculate the
    final position.  If a particle leaves its initial cell, the cell that the
    particle moves to is calculated by finding what facet the particle crosses, and
    looking up what the neighbor cell is. This is repeated until the cell
    containing the final position is found.

    \param ctrl is a DTcontrol_C object

    :cvar double dtRemaining: The time a particle has left to move in the new cell
                              the particle has entered.
    :cvar int pDim: Number of spatial coordinates in the particle location.
    :cvar double tStart: The time at which a particle starts its move in the
                                 current cell.

  */
  void move_neutral_particle_species(SegmentedArrayPair<Ptype::cartesian_xyz>& psa)
  {
    
  }
  
} // namespace dnt

#endif
