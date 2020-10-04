  //! Advance a cartesian_xy neutral-particle species by one time-increment on a mesh.
  //  template <size_t N_CELL_FACETS>
  //    void advance_neutral_species<Ptype::cartesian_xy>(py::object particle_P, py::str species_name, py::object ctrl)
  //  template <Ptype cartesian_xy, size_t N_CELL_FACETS>
  //  template <>
  template <size_t N_CELL_FACETS>
    void advance_neutral_species_cartesian_xy(py::object particle_P, py::str species_name, py::object ctrl)
  {
    using namespace pybind11::literals;
    //using PT = Ptype::cartesian_xy; not a type
        
    // std::cout << "Hello from advance_neutral_species_cartesian_xy@0" << std::endl;
    // std::cout << "This is species " << std::string(species_name) << std::endl;

    // Get attributes from the Particle_C argument
    auto pmesh_M = particle_P.attr("pmesh_M");
    auto pDim = particle_P.attr("particle_dimension").cast<int>();
    auto sap = particle_P.attr("sap_dict")[species_name].cast<SegmentedArrayPair<Ptype::cartesian_xy> *>();
    // Note that sap is a pointer.
    // Could also get sap by first casting to an std::map:
    //   auto sap_map = particle_P.attr("sap_dict").cast<std::map<std::string, SegmentedArrayPair<Ptype::cartesian_xy> *>>();
    //   auto sap = sap_map[std::string(species_name)];

    // Get call-back functions for boundary-crossing particles
    // The following lines were used when we were using the Python callbacks:
    // ?? Test for is_none() first?
    // auto bcFunctionDict = particle_P.attr("pmesh_bcs").attr("bc_function_dict").cast<py::dict>();
    // A py::dict needs a py::str key value. Convert facValue to a string first.
    // auto bcFunction = bcFunctionDict[py::str(std::to_string(facValue))].cast<py::dict>();
    // bcFunction[py::str(species_name)](ipOut, species_name, mFacet, "dx_fraction"_a=dxFraction, "facet_normal"_a=facetNormalVector);
    auto pMBCs = particle_P.attr("pmesh_bcs").cast<ParticleMeshBoundaryConditions<Ptype::cartesian_xy>>(); // (m)
    auto userPBFs = pMBCs.user_particle_boundary_functions; // Need this to call the boundary functions. // (m)
    auto bcFunctionDict = pMBCs.bc_function_dict; // (m)
    
    auto traj_T = particle_P.attr("traj_T"); // There's no cast<>() here as we only need to check if traj_T is None below.

    // Get attributes from the DTcontrol_C argument
    auto dt = ctrl.attr("dt").cast<double>(); // This does a COPY.
    auto step = ctrl.attr("timeloop_count").cast<int>();
    auto time = ctrl.attr("time").cast<double>();
    auto MAX_FACET_CROSS_COUNT = ctrl.attr("MAX_FACET_CROSS_COUNT").cast<size_t>();

    // Get particle-mesh data
    auto mesh = pmesh_M.attr("mesh").cast<dolfin::Mesh>();
    auto tDim = mesh.topology().dim();
    auto meshEntityArrays = pmesh_M.attr("mea_object").cast<MeshEntityArrays<N_CELL_FACETS> *>();
    auto particleBoundaryMarker = pmesh_M.attr("particle_boundary_marker").cast<dolfin::MeshFunction<size_t>>();
    auto NO_CELL = pmesh_M.attr("NO_CELL").cast<int>();
    auto NO_FACET = pmesh_M.attr("NO_FACET").cast<int>();
    //    std::cout << "NO_FACET " << NO_FACET << std::endl;
    
    // Get the data for the species to be advanced
    const py::ssize_t segmentLength = sap->get_segment_length();
    
    // Start a loop over the sap segments. segTuple contains (npSeg, psegIn, psegOut)
    py::tuple segTuple = sap->init_inout_loop(true); // 'true' returns the data pointers for
                                                     // the first segment pair, instead of the
                                                     // Numpy arrays.
    auto npSeg = segTuple[0].cast<py::ssize_t>();

    // Get the particle array pointers from the returned py::capsules.
    // psegIn:
    auto psegInCap = segTuple[1].cast<py::capsule>();
    //    std::cout << "The name of the capsule is " << psegInCap.name() << std::endl;
    Pstruct<Ptype::cartesian_xy>* psegIn = psegInCap; // Declares the type of psegIn
    //    std::cout << "psegIn[0].x_ = " << psegIn[0].x_ << std::endl;
    // psegOut:
    auto psegOutCap = segTuple[2].cast<py::capsule>();
    // std::cout << "The name of the capsule is " << psegOutCap.name() << std::endl;
    Pstruct<Ptype::cartesian_xy>* psegOut = psegOutCap; // Declares the type of psegOut

    // An alternative way to do this is to use the returned Numpy structured array
    // objects: call init_inout_loop() without the "true" argument above. Then use the
    // Python buffer protocol to get at the data.
    /*
    // psegIn:
    auto psegInNumpyArray = segTuple[1].cast<py::array_t<Pstruct<Ptype::cartesian_xy>,0>>();
    const auto psegInNumpyArray_info = psegInNumpyArray.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
    auto psegIn = static_cast<Pstruct<Ptype::cartesian_xy>*>(psegInNumpyArray_info.ptr); // Pointer to a particle struct in psegIn
    // psegOut:
    auto psegOutNumpyArray = segTuple[2].cast<py::array_t<Pstruct<Ptype::cartesian_xy>,0>>();
    const auto psegOutNumpyArray_info = psegOutNumpyArray.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
    auto psegOut = static_cast<Pstruct<Ptype::cartesian_xy>*>(psegOutNumpyArray_info.ptr); // Pointer to the array of particle structs.
    */
    
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
            if ((psegIn[ipIn].bitflags_ & Pstruct<Ptype::cartesian_xy>::DELETE_FLAG) != 0b0)
              {
                indexChange = true; // Particle SA indices are stale past this point
                                    // due to deletions.
                continue; // skip to the next particle in the "in" segment.
              }

            // If this "out" segment is full, get the next "out" segment, since
            // there are still some "in" particles to advance.  If there are no more
            // "out" segments, allocate a new one.
            if (ipOut == segmentLength)
              {
                // Using pointers to the plain arrays:
                py::tuple segTuple = sap->get_next_out_segment(true);                
                // auto psegOutCap = segTuple[0].cast<py::capsule>();
                // Pstruct<Ptype::cartesian_xy>* psegOut = psegInCap;
                psegOut = segTuple[0].cast<py::capsule>(); // This compresses the above 2 lines into 1
                
                // If using Numpy structured arrays:
                /*
                py::tuple segTuple = sap->get_next_out_segment();
                auto psegOutNumpyArray = segTuple[0].cast<py::array_t<Pstruct<Ptype::cartesian_xy>,0>>();
                const auto psegOutNumpyArray_info = psegOutNumpyArray.request();
                const auto psegOut = static_cast<Pstruct<Ptype::cartesian_xy>*>(psegOutNumpyArray_info.ptr);
                */
                
                ipOut = 0; // Reset the counter for the new segment
              }
            
            psegOut[ipOut] = psegIn[ipIn]; // Copy all of this particle's data from the
                                           // input slot to the output slot

            // Move the neutral particle
            // 1. Save the initial position
            psegOut[ipOut].x0_ = psegOut[ipOut].x_;
            psegOut[ipOut].y0_ = psegOut[ipOut].y_;
            // 2. Move the particle. NON-RELATIVISTIC: v and u are the same
            psegOut[ipOut].x_ += psegOut[ipOut].ux_*dt;
            psegOut[ipOut].y_ += psegOut[ipOut].uy_*dt;

            auto pCellIndex = psegOut[ipOut].cell_index_;

            //            std::cout << "advance_neutral_species. step: " << step << " ipOut: " << ipOut << " x_: " << psegOut[ipOut].x_ << std::endl;
            
            // Loop until the particle is in the current cell
            auto mLastFacet = NO_FACET;
            size_t facetCrossCount = 0;
            auto tStart = time - dt;
            auto dtRemaining = dt;

            // Create a view to this cell in the dolfin::Mesh object so that we can
            // call the entities() function to get indices of vertices, facets, etc.
            // dolfin::Cell* pcellPtr = new dolfin::Cell(mesh, pCellIndex);
            pcellPtr = new dolfin::Cell(mesh, pCellIndex);
            auto vertices = pcellPtr->entities(0);

            // std::cout << "vertices " << vertices[0] << ", " << vertices[1] << std::endl;
            
            particlePosition[0] = psegOut[ipOut].x_;
            particlePosition[1] = psegOut[ipOut].y_;
            // std::cout << "Calling is_inside_vertices" << std::endl;

            // is_inside_vertices() is declared in dolfin_functions.h.
            while (! is_inside_vertices(mesh, vertices, particlePosition))
              {
                //                std::cout << "Hello from call to is_inside_vertices" << std::endl;
                /*
                  The particle has left this cell.  We need to track it across each
                  facet in case there's a boundary-condition on that facet.
                */
                facetCrossCount += 1;
                // Check for an abnormal number of facet crossings:
                if (facetCrossCount > MAX_FACET_CROSS_COUNT)
                  {
                    std::cout << "DnT:particle.h:advance_neutral_species: !!! MAX_FACET_CROSS_COUNT exceeded!!!" << std::endl;
                    exit(EXIT_FAILURE);
                  }
                
                // Compute dx[], the move vector that starts in the current cell
                // Current particle position
                pCoord2[0] = psegOut[ipOut].x_;
                pCoord2[1] = psegOut[ipOut].y_;
                // Starting particle position in this cell
                pCoord2[3] = psegOut[ipOut].x0_;
                pCoord2[4] = psegOut[ipOut].y0_;
                // The move vector in this cell
                dx[0] = pCoord2[0] - pCoord2[3];
                dx[1] = pCoord2[1] - pCoord2[4];

                // bool returnStdArray = true;
                bool returnStdArray = false;
                py::tuple facetTupl = find_facet(pmesh_M, &pCoord2[pDim], dx, pCellIndex, returnStdArray);
                // This returns the tuple (facet, dxFraction, facet_normal_vectors[facet])
                // Extract the values from the tuple:
                auto cFacet = facetTupl[0].cast<int>();
                auto dxFraction = facetTupl[1].cast<double>();
                auto facetNormalVector = facetTupl[2].cast<py::array_t<double>>();
 
                if (cFacet != NO_FACET)
                  {
                    tStart = tStart + dxFraction*dtRemaining; // The starting time in the new cell
                    dtRemaining = (1.0 - dxFraction)*dtRemaining;

                    // Compute the crossing point
                    psegOut[ipOut].x0_ = psegOut[ipOut].x0_ + dxFraction*dx[0];
                    psegOut[ipOut].y0_ = psegOut[ipOut].y0_ + dxFraction*dx[1];
                    
                    // Look up the mesh-level index of this facet...
                    auto mFacetArray = pcellPtr->entities(tDim-1); // The mesh-level indices of the facets of this cell.
                    auto mFacet = mFacetArray[cFacet]; // The mesh-level index of the facet just crossed.

                    // mFacet should never be the same as the last facet crossed: check this
                    if (mFacet == (size_t)mLastFacet) // If the particle has crossed the same facet twice in succession, there's an error:
                      {
                        std::string errorMsg = "In advance_neutral_species(): The mesh index of the facet crossed is "  + std::to_string(mFacet) + ", the same as the last facet crossed. This should not happen since the particle cannot cross the same facet twice in one move!";
                        std::cout << errorMsg << std::endl;                        
                        exit(EXIT_FAILURE);
                      }
                    else // The particle has crossed a new facet.
                      {
                        mLastFacet = mFacet;
                      }
                    // ...and get the value of the facet marker.
                    auto facValue = particleBoundaryMarker[mFacet];
                    //                    int facValue = 0;
                    // Check if this facet has a non-zero marker, indicating that,
                    // e.g., the facet is a boundary.
                    if (facValue != 0)
                      {
                        // Call the function associated with this value.
                        if ((psegOut[ipOut].bitflags_ & Pstruct<Ptype::cartesian_xy>::TRAJECTORY_FLAG) != 0b0) // If this is a trajectory particle.
                          {
                            // Get the storage index that currently identifies this
                            // particle in the trajectory list of particles.
                            auto fullIndex = sap->get_full_index(ipIn, "in");
                            particle_P.attr("record_trajectory_datum")(std::string(species_name), ipOut, fullIndex, step, tStart, "facet_crossing"_a=true); // Does this work?
                          }
                        auto bcFunction = bcFunctionDict[std::make_pair(facValue,species_name)];
                        (userPBFs.*bcFunction)(psegOut[ipOut], species_name, mFacet, dx, dxFraction, facetNormalVector);
                      }

                    // Look up the cell index of the new cell.
                    auto pCellIndexNew = meshEntityArrays->get_cell_neighbors(pCellIndex)[cFacet];

                    /* If the particle has left the mesh, and has been deleted, end
                       the search.  If it hasn't been deleted (e.g., reflected),
                       continue tracking it.
                    */
                    if (pCellIndexNew == NO_CELL)
                      {
                        if ((psegOut[ipOut].bitflags_ & Pstruct<Ptype::cartesian_xy>::DELETE_FLAG) == 0b1)
                          {
                            psegOut[ipOut].cell_index_ = NO_CELL;
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
                else // The crossed facet is NO_FACET, which shouldn't happen.
                  {
                    std::string errorMsg = "In move_neutral_particle_species(): The cell index of the facet crossed is NO_FACET (" + std::to_string(cFacet) + "). This should not happen since the particle has left its initial cell!";
                    std::cout << errorMsg << std::endl;
                    exit(EXIT_FAILURE);
                  } // END: if (cFacet != NO_FACET)

                // The following is needed in C++ because of is_inside_vertices() usage.
                // Prepare for the next 'while' iteration: need the new cell index and the cell vertices:
                delete pcellPtr; // To avoid a leak
                pcellPtr = new dolfin::Cell(mesh, pCellIndex);
                vertices = pcellPtr->entities(0);
                
                particlePosition[0] = psegOut[ipOut].x_;
                particlePosition[1] = psegOut[ipOut].y_;
              } // END: while (! is_inside_vertices(mesh, vertices, particlePosition))
            delete pcellPtr; // Not needed past this point.

            // Record the number of facet-crossings
            psegOut[ipOut].crossings_ = facetCrossCount;

            if ((psegOut[ipOut].bitflags_ & Pstruct<Ptype::cartesian_xy>::DELETE_FLAG) == 0b0) // If this particle has not been deleted...
              {
                /*
                  If particle indices in the "out" SA have changed (i.e., indexChange
                  is 'true') due to deletions, then update this particle's index where
                  needed.  E.g., if this is a trajectory particle, update its SA
                  index in the list of trajectory particles.
                */
                if (!traj_T.is_none())
                  {
                    if ((psegOut[ipOut].bitflags_ & Pstruct<Ptype::cartesian_xy>::TRAJECTORY_FLAG) != 0b0) // If this is a trajectory particle.
                      {
                        if (indexChange == true)
                          {
                            particle_P.attr("update_trajectory_particleId")(species_name, ipIn, ipOut);
                          }
                      }
                  }
                // Advance the "out" array counter for the next particle
                ipOut += 1;
              }
            else // This particle has been deleted
              {
                indexChange = true; // This indicates that particle SA indices are changed past this point due to deletions.
                // If this was a trajectory particle, remove it's index from the trajectory-particle list.
                if (!traj_T.is_none())
                  {
                    if ((psegIn[ipIn].bitflags_ & Pstruct<Ptype::cartesian_xy>::TRAJECTORY_FLAG) != 0b0)
                      {
                        particle_P.attr("remove_trajectory_particleId")(species_name, ipIn, ipOut, step, time, dt);
                      }
                  }
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
        py::tuple segTuple = sap->get_next_segment("in", true);
        npSeg = segTuple[0].cast<py::ssize_t>();
        if (npSeg >  0)
          {
            psegIn = segTuple[1].cast<py::capsule>(); // This is a compressed expression: there's a cast performed to make psegIn, which has it's type from above (Pstruct<Ptype::cartesian_xy>* psegIn)
          }
        else
          {
            // py::print("*segTuple: ", *segTuple);
            psegIn = segTuple[1].cast<nullptr_t>(); // This will fail if segTuple[1] is not None
            assert(psegIn == nullptr);
          }

        // If we use Numpy structured arrays:
        /*
        py::tuple segTuple = sap->get_next_segment("in");
        npSeg = segTuple[0].cast<py::ssize_t>();
        auto psegInNumpyArray = segTuple[1].cast<py::array_t<Pstruct<Ptype::cartesian_xy>,0>>();
        const auto psegInNumpyArray_info = psegInNumpyArray.request();
        psegIn = static_cast<Pstruct<Ptype::cartesian_xy>*>(psegInNumpyArray_info.ptr);
        */
      } // while there are more "in" segments, keep looping

    // Set values that will be used to keep track of the particle arrays
    if (ipOut != segmentLength) // This catches the case where we exit the loop when ipOut = segmentLength and there are no more "in" segments.  Otherwise, we would add segmentLength to particleCount twice.
      {
        particleCount += ipOut;
      }
    //    std::cout << "particleCount: " << particleCount << std::endl;
    sap->set_number_of_items("out", particleCount);
  }
  // ENDDEF: void advance_neutral_species(py::object particle_P, py::str species_name, py::object ctrl)
  
