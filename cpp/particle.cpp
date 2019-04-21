/*! \file particle.h

  \brief This file has the source code for a C++ implementation of Particle storage.

  \namespace dnt

  \sa SegmentedArrayPair
*/

#include <iostream>
#include "pstruct.h"

//using namespace std; 

//template<typename PS>
//using SegList = std::vector<py::array_t<PS,0>>;

namespace dnt
{
/*! \class particle
    \brief The particle class is a C++ version of the Python class Particle_C.

    This class implements 

    \sa Segmentedarraypair

*/
  template<Ptype PT>
  class SegmentedArrayPair
  {

  private:
//    int segmentLength;
    py::ssize_t segmentLength;
  
  }
  //! Advance neutral particles by one time increment on a mesh.
  /*!

    \param psa is a SegmentedArrayPair object of the appropriated type for one
    species.

    \return void

  */
  //      using SAP = SegmentedArrayPair<PT>;
  template <typename S>
  void pstruct_to_double(S& p3D, double* point)
  {
    // The generic version
    (npSeg, psegIn, psegOut) = psa.init_inout_loop();
    
  }
  template <>
  void pstruct_to_double<DnT_pstruct1D>(DnT_pstruct1D& p1D, double* point)
  {
    point[0] = p1D.x_;
    // how about point = &(p1D.x_);
   
  }
  template <Ptype PT>
  void move_neutral_particles(SegmentedArrayPair<PT>& psa)
  {

    // Start looping over the segments
  }
  void move_neutral_particles(SegmentedArrayPair<Ptype::cartesian_x>& psa)
  {
    
  }
  
    def move_neutral_particles_CPP(self, ctrl):
        """Advance neutral particles by one time increment on a mesh.

           Compute change in position in time dt. Use an explicit method to calculate the
           final position.  If a particle leaves its initial cell, the cell that the
           particle moves to is calculated by finding what facet the particle crosses, and
           looking up what the neighbor cell is. This is repeated until the cell
           containing the final position is found.

           :param ctrl: A DTcontrol_C object

           :cvar double dtRemaining: The time a particle has left to move in the new cell
                                      the particle has entered.
           :cvar int pDim: Number of spatial coordinates in the particle location.
           :cvar double tStart: The time at which a particle starts its move in the
                                 current cell.

        """

        printInfoBoundaryCrossing = False

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        # Set local names for the passed parameters
        dt = ctrl.dt
        step = ctrl.timeloop_count
        time = ctrl.time

        pmesh_M = self.pmesh_M
        pDim = self.particle_dimension

        # Scratch space
        pCoord2 = self.pcoord2 # Can hold x,y,z, x0,y0,z0 (or a subset of these)
        dx = self.dx # dx is the distance moved in one step.

        for sn in self.neutral_species:

            # Invariant parameters
#            print 'sn = ', sn

#        nlocal = self.particles.get_species_particle_count(sp, echoFlag = False)
            
            if self.get_species_particle_count(sn) == 0: continue

            psa = self.pseg_arr[sn] # segmented array for this species

            dnt_cpp.move_neutral_particle_segment(psa)

#           Move all the particles in this species
            (npSeg, psegIn, psegOut) = psa.init_inout_loop() # (number of particles in
#            (npSeg, psegIn) = psa.init_inout_loop() # (number of particles in
                                                      # this segment, ref to
                                                      # segment)
            particleCount = 0
            # ipOut counts particles being written to the current
            # 'out' segment.
            ipOut = 0
            indexChange = False # Flag to indicate if particle SA indices have
                                # changed.  This affects, e.g., trajectories, which
                                # use SA indices to identify particles chosen for
                                # trajectory plots.
            while isinstance(psegIn, np_m.ndarray): # Keep looping until we run
                                                    # out of 'in' segments



                    
                for ipIn in range(npSeg): # Loop on the particles in this 'in'
                                           # segment
                    # psegIn[ipIn] has the 'x', 'y', 'z', 'x0',..., 'ux', 'uy',... values of ith item
                    # So psegIn[ipIn][0:3] is 'x', 'y', 'z'.
                    # Can't use slice syntax here, because the array data are not of homogeneous type.

                    # Skip deleted particles
                    if psegIn[ipIn]['bitflags'] & self.DELETE_FLAG != 0:
                        indexChange = True # Particle SA indices are stale past
                                          # this point due to deletions.
                        continue # skip to the next particle

                    # If this 'out segment is full, get the next 'out' segment,
                    # provided there are still some 'in' particles to advance.
                    if ipOut == self.SEGMENT_LENGTH:
                        psegOut = psa.get_next_out_segment()
                        ipOut = 0 # Reset the counter for the new segment
#                    if ipOut == 0: psegOut = psa.get_next_out_segment()

                    # The following COPIES data, rather than just setting a
                    # reference to existing data.
                    psegOut[ipOut] = psegIn[ipIn] # Copy this particle's data
                                                  # from the input slot to the
                                                  # output slot

# Also, check if the following copies, or resets the reference. Ans: it COPIES.
#                    p_out = psegOut[ipOut]
# To get a REFERENCE, need to use SLICING syntax:
#                    p_out = psegIn[ipIn][:]
# But you can't get a reference to a single element of psegIn[]. See:
# http://stackoverflow.com/questions/23654088/reference-of-a-single-numpy-array-element
# See also HPL p. 137. b=a[1,:], but b is still an array.
# What about getting a ref to one particle returned from a function? That's possible because it uses the stack?

                    ## Move the particle
                    # NON-RELATIVISTIC: v and u are the same
                    for coord in self.position_coordinates:
                        coord0 = coord+'0'
                        psegOut[ipOut][coord0] = psegOut[ipOut][coord] # Save the starting positions
                        ucomp = 'u'+coord
                        psegOut[ipOut][coord] += psegOut[ipOut][ucomp]*dt

#                    print "x0=", psegOut[ipOut]['x0'],"x=", psegOut[ipOut]['x']
#                    print "psegOut", psegOut[ipOut]

                    # Check if the particle has crossed into a
                    # different cell. If it has, then:
                    #   1. Find out which facet it crossed
                    #   2. Check for a boundary-condition on that facet.
                    #      a. If none (facet marker is 0), go to 3.
                    #      b. If there is one, call the BC handler. If
                    #      the BC allows this particle to continue,
                    #      goto 3. If it doesn't, the particle is
                    #      deleted.
                    #
                    #   3. Set the x0,y0,z0 coords to the crossing point.

#                    print 'ip, index =', ipOut, psegOut[ipOut]['cell_index']
                    pCellIndex = psegOut[ipOut]['cell_index']
#                    print fncName, ": ip, pindex", ipOut, pCellIndex, "cell index:", pmesh_M.compute_cell_index(psegOut[ipOut])

                    # Loop until the particle is in the current cell
                    mLastFacet = pmesh_M.NO_FACET
                    facetCrossCount = 0
                    tStart = time - dt
                    dtRemaining = dt
#                    while not pmesh_M.is_inside(psegOut[ipOut], pCellIndex):
                    while not pmesh_M.is_inside_CPP(psegOut[ipOut], pCellIndex):
                        # The particle has left this cell.  We
                        # need to track it across each facet in case
                        # there's a boundary-condition on that facet.
#                        print fncName, "particle has migrated"

                        facetCrossCount += 1
                        # Check for an abnormal number of facet crossings:
                        if facetCrossCount > self.MAX_FACET_CROSS_COUNT:
                            errorMsg = "%s !!! MAX_FACET_CROSS_COUNT exceeded!!!" % (fncName)
                            sys.exit(errorMsg)

                        # Compute dx[], the move vector that starts in
                        # the current cell
                        i=0
                        for coord in self.position_coordinates:
                            pCoord2[i] = psegOut[ipOut][coord] # Present position
                            coord0 = coord+'0'
                            pCoord2[i+pDim] = psegOut[ipOut][coord0] # Start of path in this cell
                            i+=1
#                        pCoord2[:] = psegOut[ipOut] # Alias to the position coordinates
                        dx = pCoord2[0:pDim] - pCoord2[pDim:2*pDim] # Move vector

# could return the crossing-point in pCoord2[]
                        (cFacet, dxFraction, facetNormal) = pmesh_M.find_facet(pCoord2[pDim:2*pDim], dx, pCellIndex)
                        if cFacet != pmesh_M.NO_FACET:
#                            print "facet crossed is", cFacet
                            tStart = tStart + dxFraction*dtRemaining # The starting time in the new cell
                            dtRemaining = (1.0 - dxFraction)*dtRemaining

                            # Compute the crossing point
#                           print fncName, "Before truncation p =", psegOut[ipOut]
                            i=0 # i indexes the move-vector coordinates (dx[i])
                            for coord in self.position_coordinates:
                                coord0 = coord+'0'
                                psegOut[ipOut][coord0] = psegOut[ipOut][coord0] + dxFraction*dx[i]
                                i+=1
#                            print "After truncation p =", psegOut[ipOut]

                            # Look up the mesh-level index of this facet...
                            mFacet = pmesh_M.cell_entity_index_dict['facet'][pCellIndex][cFacet]
# mFacet should never be the same as the last facet crossed: check this
                            if mFacet == mLastFacet: # If the particle has crossed the same facet twice in succession...
                                errorMsg = "%s The mesh index of the facet crossed is %d, the same as the last facet crossed. This should not happen since the particle cannot cross the same facet twice in one move!" % (fncName, mFacet)
                                sys.exit(errorMsg)
                            else: # The particle has crossed a new facet.
                                mLastFacet = mFacet
                            # ...and get the value of the facet marker.
# mFacet is a numpy.uint32 (size_t), but the FacetFunction wants an int argument.
#                            print "type is:", type(mFacet)
                            facValue = pmesh_M.particle_boundary_marker[mFacet]
                            # Check if this facet has a non-zero marker, indicating that,
                            # e.g., the facet is a boundary.
                            if facValue != 0:
                                # Call the function associated with this value.

                                if psegOut[ipOut]['bitflags'] & self.TRAJECTORY_FLAG != 0: # If this is a trajectory particle.

                                    if printInfoBoundaryCrossing is True: print("Recording boundary-crossing for particle", psegOut[ipOut]['unique_ID'])

                                    # Get the storage index that currently identifies this
                                    # particle in the trajectory list of particles.
                                    fullIndex = psa.get_full_index(ipIn, 'in')
                                    self.record_trajectory_datum(sn, psegOut[ipOut], fullIndex, step, tStart, facet_crossing=True)
                                # A reference to dx[] is available in the BC function class.
                                self.pmesh_bcs.bc_function_dict[facValue][sn](psegOut[ipOut], sn, mFacet, dx_fraction=dxFraction, facet_normal=facetNormal)
                            # Look up the cell index of the new cell.
                            pCellIndexNew = pmesh_M.cell_neighbor_dict[pCellIndex][cFacet]

                            # If the particle has left the mesh, and has been deleted, end
                            # the search.  If it hasn't been deleted (e.g., reflected),
                            # continue tracking it.
                            if pCellIndexNew == pmesh_M.NO_CELL:
                                if psegOut[ipOut]['bitflags'] & self.DELETE_FLAG == 1:
                                    psegOut[ipOut]['cell_index'] = pmesh_M.NO_CELL
                                    break # Breaks out of the facet-crossing 'while' loop.
                                # else: The boundary did not absorb the particle, so the
                                # particle is now at the cell boundary and the cell index
                                # is unchanged.
                            else:
                                psegOut[ipOut]['cell_index'] = pCellIndexNew
                                pCellIndex = pCellIndexNew # Needed for the next iteration of the while loop.

                        else: # The crossed faced is NO_FACET, which shouldn't happen.
                            errorMsg = "%s The cell index of the facet crossed is %d. This should not happen since the particle has left its initial cell cell!" % (fncName, cFacet)
                            sys.exit(errorMsg)
#                       END:if cFacet != pmesh_M.NO_FACET:
#                   END:while not pmesh_M.is_inside(psegOut[ipOut], pCellIndex)

                    # Record the number of facet-crossings
                    psegOut[ipOut]['crossings'] = facetCrossCount
                    
                    # Check that this particle has not been deleted before
                    # incrementing the 'out' particle counter.
                    if psegOut[ipOut]['bitflags'] & self.DELETE_FLAG == 0: # If this particle has not been deleted...
                        # If particle indices in the 'out' SA have changed (i.e.,
                        # indexChange is True) due to deletions, then update this
                        # particle's index where needed.  E.g., if this is a
                        # trajectory particle, update its SA index in the list of
                        # trajectory particles.
                        if self.traj_T is not None:
                            if psegOut[ipOut]['bitflags'] & self.TRAJECTORY_FLAG != 0: # If this is a trajectory particle.
                                if indexChange is True:
#                                    print "mover: indexChange is true: updating the trajectory particle ID for particle", psegIn[ipIn]['unique_ID'], "ipIn =", ipIn, "ipOut =", ipOut
                                    self.update_trajectory_particleId(sn, ipIn, ipOut)

                        # Advance the 'out' array counter for the next particle
                        ipOut += 1
                    else: # This particle has been deleted
                        indexChange = True # This indicates that particle SA indices
                                           # are changed past this point due to
                                           # deletions.
                        # If this was a trajectory particle, remove it's index from
                        # the trajectory-particle list.
                        if self.traj_T is not None:
                            if psegIn[ipIn]['bitflags'] & self.TRAJECTORY_FLAG != 0:
#                                print "mover: Removing particle with ID", psegIn[ipIn]['unique_ID'], "ipIn", ipIn, "from trajectories", ", ipOut is", ipOut
                                self.remove_trajectory_particleId(sn, ipIn, psegOut[ipOut], step, time, dt)

                    # Check if we've reached the end of this segment.  If
                    # so, we need to start writing on a new segment.
                    if (ipOut == self.SEGMENT_LENGTH):
                        particleCount += self.SEGMENT_LENGTH
#                        ipOut = 0 # This will cause get_next_out_segment() to be called
                                   # above, if there are more 'in' particles to be processed.

                # Done with this segment.
                # Get the next one, if it exists.
                (npSeg, psegIn) = psa.get_next_segment('in')
            # End of while isinstance(psegIn, np_m.ndarray)                                   
                
            # Set values that will be used to keep track of the particle arrays
            if (ipOut != self.SEGMENT_LENGTH): # This catches the case where we exit the loop
                                               # when ipOut = SEGMENT_LENGTH and there are no
                                               # more 'in' segments.  Otherwise, we would add
                                               # SEGMENT_LENGTH to particleCount twice.
                particleCount += ipOut
            psa.set_number_of_items('out', particleCount)
            # Loop over segmented array ends

# Compute new density here?

        # Loop over Species ends
        return
#    def move_neutral_particles_CPP(self, ctrl):ENDDEF
