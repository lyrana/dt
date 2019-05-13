/*! \file particle.h

  \brief This file has the source code for a C++ implementation of particle advance.

  \namespace dnt

  \sa SegmentedArrayPair
*/

#include <iostream>
//#include "pstruct.h"
#include "dolfin.h"

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
  //  void move_neutral_particles_species(SegmentedArrayPair<PT>& psa, dolfin::Mesh& mesh, py::dict cell_vertex_dict)
  void move_neutral_particle_species(py::object ctrl, SegmentedArrayPair<PT>& psa, dolfin::Mesh& mesh, py::dict cell_vertex_dict, py::object particle_P)
  {

    //    HERE
      
    double dt = ctrl.attr("dt").cast<double>(); // Don't know if this works yet. Obviously, this does a COPY.

    int MAX_FACET_CROSS_COUNT = ctrl.attr("MAX_FACET_CROSS_COUNT").cast<int>()
    
    auto mesh_coordinates = mesh.coordinates();
    
    // Start a loop over the segments
    py::tuple stupl = psa.init_inout_loop(true); // 'true' returns the data pointers for the
                                                 // first segment pair.
    auto npSeg = stupl[0].cast<py::ssize_t>();
    // For Numpy array return: psegInArray = stupl[1].cast<py::array_t<Pstruct<PT>,0>>();
    auto psegIn = stupl[1].cast<Pstruct<PT>*>(); // the segment being read from.
    // For Numpy array return: psegOutArray = stupl[2].cast<py::array_t<Pstruct<PT>,0>>();
    auto psegOut = stupl[2].cast<Pstruct<PT>*>(); // the segment being written to.
    
    py::ssize_t particleCount = 0; // Counts the number of still-active particles for this
                                   // species.
    
    // ipOut counts particles being written to the current 'out' segment.
    py::ssize_t ipOut = 0;
    bool indexChange = false; // Flag to indicate if the particle SA indices have
                              // changed. This affects, e.g., trajectories, which use SA
                              // indices to identify particles chosen for trajectory
                              // plots.
    
    while (psegIn != nullptr) // Keep looping until we run out of 'in' segments.
      {
        for (auto ipIn = 0; ipIn < npSeg; ipIn++)
          {
            // Skip deleted particles
            //            if (psegIn[ipIn]['bitflags'] & self.DELETE_FLAG != 0)
            if (psegIn[ipIn].bitflags_ & Pstruct<PT>.DELETE_FLAG != 0b0)
              {
                indexChange = true; // Particle SA indices are stale past this point
                                    // due to deletions.
                continue; // skip to the next particle
              }

            // If this 'out' segment is full, get the next 'out' segment, since
            // there are still some 'in' particles to advance.  If there are no more
            // 'out' segments, allocate a new one.
            if (ipOut == psa.get_segment_length())
              {
                psegOut = psa.get_next_out_segment(true); // 'true' returns a data pointer
                ipOut = 0; // Reset the counter for the new segment
              }
            
            psegOut[ipOut] = psegIn[ipIn]; // Copy all of this particle's data from the
                                           // input slot to the output slot

            // Move the particle
            // NON-RELATIVISTIC: v and u are the same
            psegOut[ipOut].x0_= psegOut[ipOut].x_;  // Save the starting positions
            psegOut[ipOut].x_ += psegOut[ipOut].ux_*dt;

            /*            
              Check if the particle has crossed into a different cell. If it has,
              then:
              1. Find out which facet it crossed
              2. Check for a boundary-condition on that facet.
                 a. If none (facet marker is 0), go to 3.
                 b. If there is one, call the BC handler. If the BC allows this
                    particle to continue, goto 3. If it doesn't, the particle is
                    deleted.
              3. Set the x0,y0,z0 coords to the crossing point.
            */

            int pCellIndex = psegOut[ipOut].cell_index_;

            // Loop until the particle stays in the current cell
            int mLastFacet = Mesh_C_NO_FACET;
            int facetCrossCount = 0;
            double tStart = time - dt;
            double dtRemaining = dt;

            
            // while (! pmesh_M.is_inside_CPP(psegOut[ipOut], pCellIndex))
            //            py::array_t<int> vertices = cell_vertex_dict[pCellIndex];
            //            while (! cell_contains_point(mesh, cell_vertex_dict, psegOut[ipOut]))
            while (! particle_stays_in_cell(psegOut[ipOut], mesh, cell_vertex_dict))
              {
                /*
                  The particle has left this cell.  We need to track it across each
                  facet in case there's a boundary-condition on that facet.
                */
                facetCrossCount += 1;
                // Check for an abnormal number of facet crossings:
                if (facetCrossCount > MAX_FACET_CROSS_COUNT)
                  {
                    errorMsg = "%s !!! MAX_FACET_CROSS_COUNT exceeded!!!" % (fncName);
                    sys.exit(errorMsg);
                  }

                // Compute dx[], the move vector that starts in the current cell
                int i=0;
                for (coord in self.position_coordinates)
                  {
                    pCoord2[i] = psegOut[ipOut][coord];// Present position
                    coord0 = coord+'0';
                    pCoord2[i+pDim] = psegOut[ipOut][coord0]; // Start of path in this cell
                    i+=1;
                  }
                dx = pCoord2[0:pDim] - pCoord2[pDim:2*pDim]; // Move vector

                (cFacet, dxFraction, facetNormal) = pmesh_M.find_facet(pCoord2[pDim:2*pDim], dx, pCellIndex);
                if (cFacet != Mesh_C_NO_FACET)
                  {
                    tStart = tStart + dxFraction*dtRemaining; // The starting time in the new cell
                    dtRemaining = (1.0 - dxFraction)*dtRemaining;

                    // Compute the crossing point
                    int i=0; // i indexes the move-vector coordinates (dx[i])
                    for (coord in self.position_coordinates)
                      {
                        coord0 = coord+'0';
                        psegOut[ipOut][coord0] = psegOut[ipOut][coord0] + dxFraction*dx[i];
                        i+=1;
                      }
                    // Look up the mesh-level index of this facet...
                    mFacet = pmesh_M.cell_entity_index_dict['facet'][pCellIndex][cFacet]; // mFacet should never be the same as the last facet crossed: check this
                    if (mFacet == mLastFacet) // If the particle has crossed the same facet twice in succession...
                      {
                        errorMsg = "%s The mesh index of the facet crossed is %d, the same as the last facet crossed. This should not happen since the particle cannot cross the same facet twice in one move!" % (fncName, mFacet);
                        sys.exit(errorMsg);
                      }
                    else // The particle has crossed a new facet.
                      {
                        mLastFacet = mFacet;
                      }
                    // ...and get the value of the facet marker.
                    int facValue = pmesh_M.particle_boundary_marker[mFacet];
                    // Check if this facet has a non-zero marker, indicating that,
                    // e.g., the facet is a boundary.
                    if (facValue != 0)
                      {
                        // Call the function associated with this value.

                        if (psegOut[ipOut].bitflags_ & Pstruct<PT>.TRAJECTORY_FLAG != 0b0) // If this is a trajectory particle.
                          {
                            // Get the storage index that currently identifies this
                            // particle in the trajectory list of particles.
                            fullIndex = psa.get_full_index(ipIn, 'in');
                            self.record_trajectory_datum(sn, psegOut[ipOut], fullIndex, step, tStart, facet_crossing=True);
                          }
                        // A reference to dx[] is available in the BC function class.
                        self.pmesh_bcs.bc_function_dict[facValue][sn](psegOut[ipOut], sn, mFacet, dx_fraction=dxFraction, facet_normal=facetNormal);
                      }
                    // Look up the cell index of the new cell.
                    pCellIndexNew = pmesh_M.cell_neighbor_dict[pCellIndex][cFacet];

                    /* If the particle has left the mesh, and has been deleted, end
                       the search.  If it hasn't been deleted (e.g., reflected),
                       continue tracking it.
                    */
                    if (pCellIndexNew == pmesh_M.NO_CELL)
                      {
                        if (psegOut[ipOut].bitflags_ & Pstruct<PT>.DELETE_FLAG == 1)
                          {
                            psegOut[ipOut].cell_index_ = pmesh_M.NO_CELL
                              break; // Breaks out of the facet-crossing 'while' loop.
                          }
                        // else: The boundary did not absorb the particle, so the
                        // particle is now at the cell boundary and the cell index
                        // is unchanged.
                      }
                    else
                      {
                        psegOut[ipOut].cell_index_ = pCellIndexNew;
                        pCellIndex = pCellIndexNew; // Needed for the next iteration of the while loop.
                      }
                  }
                else // The crossed faced is NO_FACET, which shouldn't happen.
                  {
                    std::string errorMsg = "In move_neutral_particle_species(): The cell index of the facet crossed is" + std::to_string(cFacet) + ". This should not happen since the particle has left its initial cell cell!";
                    dnt_error(errorMsg);
                  } // END: if (cFacet != Mesh_C_NO_FACET)
              } // END: while (! is_inside())

            // Record the number of facet-crossings
            psegOut[ipOut].crossings_ = facetCrossCount;
                    
            // Check that this particle has not been deleted before incrementing the
            // 'out' particle counter.
            if (psegOut[ipOut].bitflags_ & Pstruct<PT>.DELETE_FLAG == 0b0) // If this particle has not been deleted...
              {
                /*
                  If particle indices in the 'out' SA have changed (i.e., indexChange
                  is 'true') due to deletions, then update this particle's index where
                  needed.  E.g., if this is a trajectory particle, update its SA
                  index in the list of trajectory particles.
                */
                if (traj_T != nullptr)
                  {
                    if psegOut[ipOut].bitflags_ & Pstruct<PT>.TRAJECTORY_FLAG != 0b0: // If this is a trajectory particle.
                    if (indexChange is True)
                      {
                        self.update_trajectory_particleId(sn, ipIn, ipOut);
                      }
                  }
                // Advance the 'out' array counter for the next particle
                ipOut += 1;
              }
            else // This particle has been deleted
              {
                indexChange = true; // This indicates that particle SA indices are changed past this point due to deletions.
                // If this was a trajectory particle, remove it's index from the trajectory-particle list.
                if (self.traj_T is not None)
                  {
                    if (psegIn[ipIn].bitflags_ & Pstruct<PT>.TRAJECTORY_FLAG != 0b0)
                      {
                        self.remove_trajectory_particleId(sn, ipIn, psegOut[ipOut], step, time, dt);
                      }
                  }
              }
            /*
              Check if we've reached the end of this segment.  If so, we need to start
              writing on a new segment.  If there are no more segments, allocate a new one.
            */
            if (ipOut == self.SEGMENT_LENGTH)
              {
                particleCount += self.SEGMENT_LENGTH;
              }

          } // loop over particles in the 'in' segment
        // Done with this segment.Get the next one, if it exists.
        (npSeg, psegIn) = psa.get_next_segment('in');
      } // while there are more 'in' segments

    // Set values that will be used to keep track of the particle arrays
    if (ipOut != self.SEGMENT_LENGTH) // This catches the case where we exit the loop when ipOut = SEGMENT_LENGTH and there are no more 'in' segments.  Otherwise, we would add SEGMENT_LENGTH to particleCount twice.
      {
        particleCount += ipOut;
      }
    psa.set_number_of_items('out', particleCount);
  }
  // void move_neutral_particles(SegmentedArrayPair<PT>& psa) ENDDEF



  
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
  void move_neutral_particle_species(SegmentedArrayPair<Ptype::cartesian_x>& psa)
  {
    
  }
  
