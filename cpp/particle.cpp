/*! \file particle.h

  \brief This file has the source code for a C++ implementation of particle advance.

  \namespace dnt

  \sa SegmentedArrayPair
*/

#include <iostream>
//#include "pstruct.h"
#include "SegmentedArrayPair.h"

//#include "dolfin.h"

//using namespace std; 

//template<typename PS>
//using SegList = std::vector<py::array_t<PS,0>>;

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

    // Get parameters needed from ctrl arg
    auto dt = ctrl.attr("dt").cast<double>(); // This does a COPY.
    auto MAX_FACET_CROSS_COUNT = ctrl.attr("MAX_FACET_CROSS_COUNT").cast<int>();

    // Get attributes needed from Particle_C arg
    auto pmesh_M = particle_P.attr("pmesh_M");
    auto pDim = particle_P.attr("particle_dimension").cast<int>();
    auto pseg_arr = particle_P.attr("pseg_arr").cast<py::dict>();
    
    // Get the data for the species to be advanced
    //    auto pseg_arr_map = particle_P.attr("pseg_arr").cast<std::map<std::string, SegmentedArrayPair<PT *>>();
    //    auto psa = pseg_arr_map[std::string(species_name)];

    // Access the SAP without constructing a std::map. Note that psa is a pointer.
    auto psa = pseg_arr[species_name].cast<SegmentedArrayPair<PT> *>();

    const py::ssize_t segmentLength = psa->get_segment_length();
    //    auto mesh_coordinates = mesh.coordinates();
    
    // Start a loop over the psa segments
    py::tuple stupl = psa->init_inout_loop(true); // 'true' returns the data pointers for
                                                 // the first segment pair, instead of the
                                                 // Numpy arrays.
    auto npSeg = stupl[0].cast<py::ssize_t>();

    // For Numpy array return: psegInArray = stupl[1].cast<py::array_t<Pstruct<PT>,0>>();
    auto psegIn = stupl[1].cast<Pstruct<PT>*>(); // the segment being read from.

    // For Numpy array return: psegOutArray = stupl[2].cast<py::array_t<Pstruct<PT>,0>>();
    auto psegOut = stupl[2].cast<Pstruct<PT>*>(); // the segment being written to.
    
    py::ssize_t particleCount = 0; // Counts the number of still-active particles for this
                                   // species.
    
    // ipOut counts particles being written to the current "out" segment.
    py::ssize_t ipOut = 0;
    bool indexChange = false; // Flag to indicate if the particle SA indices have
                              // changed. This affects, e.g., trajectories, which use SA
                              // indices to identify particles chosen for trajectory
                              // plots.
    
    while (psegIn != nullptr) // Keep looping until we run out of "in" segments.
      {
        for (auto ipIn = 0; ipIn < npSeg; ipIn++)
          {
            // Skip deleted particles
            //            if (psegIn[ipIn]['bitflags'] & self.DELETE_FLAG != 0)
            if (psegIn[ipIn].bitflags_ & Pstruct<PT>::DELETE_FLAG != 0b0)
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
                psegOut = psa->get_next_out_segment(true); // 'true' returns a data pointer
                                                          // instead of a Numpy array.
                ipOut = 0; // Reset the counter for the new segment
              }
            
            psegOut[ipOut] = psegIn[ipIn]; // Copy all of this particle's data from the
                                           // input slot to the output slot

            /*
              Check if we've reached the end of this segment.  If so, we need to start
              writing on a new segment.  If there are no more segments, allocate a new one.
            */
            if (ipOut == segmentLength)
              {
                particleCount += segmentLength;
              }

          } // loop over particles in the "in" segment
        // Done with this segment.Get the next one, if it exists.
        (npSeg, psegIn) = psa->get_next_segment("in");
      } // while there are more "in" segments, keep looping

    // Set values that will be used to keep track of the particle arrays
    if (ipOut != segmentLength) // This catches the case where we exit the loop when ipOut = segmentLength and there are no more "in" segments.  Otherwise, we would add segmentLength to particleCount twice.
      {
        particleCount += ipOut;
      }
    psa->set_number_of_items("out", particleCount);
  }
  // ENDDEF: void move_neutral_particles(SegmentedArrayPair<PT>& psa)



  
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
  void move_neutral_particle_species(SegmentedArrayPair<Ptype::cartesian_x_y_z>& psa)
  {
    
  }
  
} // namespace dnt
