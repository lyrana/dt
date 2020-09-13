  //! Advance a spherical_r charged-particle species by one time-increment on a mesh.
  /*!          

    Apply interpolated electric force to charged particles.  Compute the change in
    velocity and position in time dt. Use an explicit method to integrate the orbit.

    If a particle leaves its initial cell, the cell that the particle moves to is
    calculated by finding what facet the particle crosses, and looking up what the
    neighbor cell is. This is repeated until the cell containing the final position
    is found.

    \param[in,out] particle_P is a ref to a Particle_C Python object.
    \param[in] species_name is the py::str name of the species to be advanced.
    \param[in] ctrl is a ref to a DTcontrol_C Python object.

    \var SegmentedArrayPair<PT>* sap
    This is a pointer-to-SegmentedArrayPair containing the particle data
    for one neutral species.
    \var double dtRemaining
    The time a particle has left to move in the new cell the particle has entered.
    \var int pDim
    Number of spatial coordinates in the particle location.
    \var double tStart
    The time at which a particle starts its move in the current cell.

  */
  // template <Ptype PT, size_t TOPOL_DIM>
  //self.integrators[sn](self, sn, ctrl, neg_E_field=neg_E_field, external_E_field=external_E_field, accel_only=accel_only)
  
  //  template <Ptype PT, size_t N_CELL_FACETS>
    //    void advance_charged_species_in_E_field(py::object particle_P, py::str species_name, py::object ctrl, dolfin::Function& neg_E_field, dolfin::Function& external_E_field, bool accel_only)
//  template <Ptype PT, size_t N_CELL_FACETS>
//    void advance_charged_species_in_E_field(py::object particle_P, py::str species_name, py::object ctrl, dolfin::Function* neg_E_field, dolfin::Function* external_E_field, bool accel_only)
//  {}

  template <size_t N_CELL_FACETS>
    void advance_charged_species_in_E_field_spherical_r(py::object &particle_P, py::str species_name, py::object &ctrl, dolfin::Function *neg_E_field, dolfin::Function *external_E_field, bool accel_only)
  {
    using namespace pybind11::literals;

    //std::cout << "Hello from advance_charged_species_in_E_field_spherical_r@1" << std::endl;
    //std::cout << "This is species " << std::string(species_name) << std::endl;

    // Maybe all of the stuff that doesn't change from call to call should be moved to "initialize", and used to initialize static variables?
    // Marked these with (m) below
    // For parallel runs, need to call "intialize" again after a rebalance.
    
    // Get attributes from the Particle_C argument
    
    auto pmesh_M = particle_P.attr("pmesh_M"); // (m)
    auto pDim = particle_P.attr("particle_dimension").cast<int>(); // (m)
    auto sap = particle_P.attr("sap_dict")[species_name].cast<SegmentedArrayPair<Ptype::spherical_r> *>();
    // Note that sap is a pointer.
    // Could also get sap by first casting to a std::map:
    //   auto sap_map = particle_P.attr("sap_dict").cast<std::map<std::string, SegmentedArrayPair<PT> *>>();
    //   auto sap = sap_map[std::string(species_name)];
    //   auto sapMap = particle_P.attr("sap_dict").cast<std::map<std::string, SegmentedArrayPair<Ptype::spherical_r> *>>();
    auto qom = particle_P.attr("qom")[species_name].cast<double>();
    // Get call-back functions for boundary-crossing particles
    // The following line was used when we were using the Python callbacks:
    // auto bcFunctionDict = particle_P.attr("pmesh_bcs").attr("bc_function_dict").cast<std::map<int, py::dict>>();
    auto pMBCs = particle_P.attr("pmesh_bcs").cast<ParticleMeshBoundaryConditions<Ptype::spherical_r>>(); // (m)
    auto userPBFs = pMBCs.user_particle_boundary_functions; // Need this to call the boundary functions. // (m)
    auto bcFunctionDict = pMBCs.bc_function_dict; // (m)
    
    // These are scratch Numpy arrays for field interpolation to the particles. They have
    // "shapes" that depend on how they were declared.
    auto negE = particle_P.attr("negE").cast<py::array_t<double>>(); // (m)
    auto Eext = particle_P.attr("Eext").cast<py::array_t<double>>(); // (m)
    auto zeroE = particle_P.attr("zeroE").cast<py::array_t<double>>(); // (m)
    // These are used for holding the fields of trajectory particles
    auto E1 = particle_P.attr("E1_comps").cast<py::array_t<double>>(); // (m)
    auto Eext1 = particle_P.attr("Eext1_comps").cast<py::array_t<double>>(); // (m)

    // To use the buffer protocol instead of a proxy
    /* py::buffer_info negEinfo = negE.request(); // request() returns metadata about the array (ptr, ndim, size, shape) */
    /* negE_array = (double*) negEinfo.ptr; // A simple C++ array */
    /* E_array_len = negEinfo.shape[0]; */
    /* E_array_ncomps = negEinfo.shape[1]; */
    // Could make an Eigen matrix from this info.
    
    auto traj_T = particle_P.attr("traj_T"); // There's no cast<>() here as we only need to check if traj_T is None below.

    // Get particle-mesh data
    auto mesh = pmesh_M.attr("mesh").cast<dolfin::Mesh>(); // (m)
    auto tDim = mesh.topology().dim(); // (m)
    auto meshEntityArrays = pmesh_M.attr("mea_object").cast<MeshEntityArrays<N_CELL_FACETS> *>(); // (m)
    auto particleBoundaryMarker = pmesh_M.attr("particle_boundary_marker").cast<dolfin::MeshFunction<size_t>>(); // (m)
    auto NO_CELL = pmesh_M.attr("NO_CELL").cast<int>(); // (m)
    auto NO_FACET = pmesh_M.attr("NO_FACET").cast<int>(); // (m)
    // std::cout << "NO_FACET " << NO_FACET << std::endl;
    
    // Get attributes from the DTcontrol_C argument
    auto dt = ctrl.attr("dt").cast<double>(); // This does a COPY.
    auto step = ctrl.attr("timeloop_count").cast<int>();
    auto time = ctrl.attr("time").cast<double>();
    auto MAX_FACET_CROSS_COUNT = ctrl.attr("MAX_FACET_CROSS_COUNT").cast<size_t>(); // (m)
    auto applySolvedElectricField = ctrl.attr("apply_solved_electric_field"); // (m)
    
    
    // Check if self-electric field should be applied to this species
    bool applySolvedElectricFieldFlag;
    std::map<std::string, bool> applySolvedElectricFieldMap;
    if (applySolvedElectricField.is_none())
      {
        applySolvedElectricFieldFlag = true;
      }
    else
      {
        applySolvedElectricFieldFlag = false;
        applySolvedElectricFieldMap = applySolvedElectricField.cast<std::map<std::string, bool>>();
        //std::cout << "advance_charged_species_in_E_field_spherical_r: applySolvedElectricFieldMap[species_name] = " << applySolvedElectricFieldMap[species_name] << std::endl;
      }

    auto applyRandomExternalElectricField = ctrl.attr("apply_random_external_electric_field");
    std::map<std::string, bool> applyRandomExternalElectricFieldMap;
    // Check if the solved electric field should be applied to this species
    bool applyRandomExternalElectricFieldFlag;
    if (applyRandomExternalElectricField.is_none())
      {
        applyRandomExternalElectricFieldFlag = true;
      }
    else
      {
        applyRandomExternalElectricFieldFlag = false;
        applyRandomExternalElectricFieldMap = applyRandomExternalElectricField.cast<std::map<std::string, bool>>();
        //std::cout << "advance_charged_species_in_E_field_spherical_r: applyRandomExternalElectricFieldMap[species_name] = " << applyRandomExternalElectricFieldMap[species_name] << std::endl;
      }
    
    // Invariant parameters used in the particle-advance loop
    auto qmdt = qom*dt;

    // Get the particles belonging to this species
    
    const py::ssize_t segmentLength = sap->get_segment_length();
    // Start a loop over the sap segments. segTuple contains (npSeg, psegIn, psegOut)
    py::tuple segTuple = sap->init_inout_loop(true); // 'true' returns the data pointers for
                                                     // the first segment pair, instead of the
                                                     // Numpy arrays.
    auto npSeg = segTuple[0].cast<py::ssize_t>(); // The number of active particles in the segment.

    // Get the particle array pointers from the py::capsules in segTuple
    // psegIn:
    auto psegInCap = segTuple[1].cast<py::capsule>();
    //    std::cout << "The name of the capsule is " << psegInCap.name() << std::endl;
    Pstruct<Ptype::spherical_r>* psegIn = psegInCap; // Declares the type of psegIn
    // std::cout << "psegIn[0].r_ = " << psegIn[0].r_ << std::endl;
    // psegOut:
    auto psegOutCap = segTuple[2].cast<py::capsule>();
    // std::cout << "The name of the capsule is " << psegOutCap.name() << std::endl;
    Pstruct<Ptype::spherical_r>* psegOut = psegOutCap; // Declares the type of psegOut

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
    auto psegOut = static_cast<Pstruct<PT>*>(psegOutNumpyArray_info.ptr); // Pointer to the array of particle structs.
    */
    
    py::ssize_t particleCount = 0; // Counts the number of still-active particles for this
                                   // species.
    
    // ipOut counts particles being written to the current "out" segment.
    py::ssize_t ipOut = 0;
    bool indexChange = false; // Flag to indicate if the particle SA indices have
                              // changed. This affects, e.g., trajectories, which use SA
                              // indices to identify particles chosen for trajectory
                              // plots.

    //std::cout << "Hello from advance_charged_species_in_E_field_spherical_r@2" << std::endl;
    
    py::array_t<double>* EsegPtr = nullptr;
    while (psegIn != nullptr) // Keep looping until we run out of "in" segments.
      {
        //std::cout << "Hello from advance_charged_species_in_E_field_spherical_r@3" << std::endl;
        // Compute the electric field at each particle

        // Eseg = None # This is for the case where no fields are to be applied to the particles,
        //        if (!neg_E_field.is_none())
        if (neg_E_field != nullptr)
          {
            //std::cout << "Hello from advance_charged_species_in_E_field_spherical_r@3.1" << " applySolvedElectricFieldFlag = " << applySolvedElectricFieldFlag << std::endl;
            if ((applySolvedElectricFieldFlag == true) || (applySolvedElectricFieldMap[species_name] == true))
              //if (applySolvedElectricFieldFlag == true)
              {
                //std::cout << "Hello from advance_charged_species_in_E_field_spherical_r@3.2" << " npSeg = " << npSeg << std::endl;
                // # std::cout << "interpolate_field_to_points" << std::endl;
                interpolate_field_to_points<Ptype::spherical_r>(neg_E_field, psegIn, npSeg, negE); // From dolfin_functions.h
                //std::cout << "Hello from advance_charged_species_in_E_field_spherical_r@4" << std::endl;
                EsegPtr = &negE;
                auto EsegProxy = EsegPtr->mutable_unchecked<2>();
                auto nComps = EsegProxy.shape(1);
                // Flip the sign get the E-field (or do this in the acceleration loop below.)
                for (ssize_t i = 0; i < npSeg; i++)
                  {
                    for (ssize_t j = 0; j < nComps; j++)
                      {
                        EsegProxy(i,j) *= -1.0;
                        //std::cout << "Eseg i, j " << i << " , " << j << " = " <<  EsegProxy(i,j) << std::endl;
                      }
                  }
              }
          }
        else
          {
            //std::cout << "Hello from advance_charged_species_in_E_field_spherical_r@5" << std::endl;
            EsegPtr = &zeroE;
          }
        // Add an external E field, if present, to E.
        if (external_E_field != nullptr)
          {
            if ((applyRandomExternalElectricFieldFlag == true) || (applyRandomExternalElectricFieldMap[species_name] == true))
              {
                // Interpolate the external electric field to the particle positions
                std::cout << "advance_charged_species_in_E_field_spherical_r: This function needs to be written" << std::endl;                
                // interpolate_random_field_to_points<Ptype::spherical_r>(external_E_field, psegIn, npSeg, Eext);
                auto EsegProxy = EsegPtr->mutable_unchecked<2>();
                auto EextProxy = Eext.mutable_unchecked<2>();
                auto nComps = EextProxy.shape(1);
                // Add the interpolated external field to the interpolated solved electric field.
                for (ssize_t i = 0; i < npSeg; i++)
                  {
                    for (ssize_t j = 0; j < nComps; j++)
                      {
                        EsegProxy(i,j) += EextProxy(i,j);
                        //std::cout << "Eseg particle, comp: " << i << " , " << j << " = " <<  EsegProxy(i,j) << std::endl;
                      }
                  }
              }
          }

        // Advance the particles in this segment
        for (auto ipIn = 0; ipIn < npSeg; ipIn++)
          {
            // If this "out" segment is full, get the next "out" segment, since
            // there are still some "in" particles to advance.  If there are no more
            // "out" segments, allocate a new one.
            // # std::cout << "if ipOut == segmentLength" << std::endl;
            if (ipOut == segmentLength)
              {
                // Using pointers to the plain arrays:
                // # std::cout << "segTuple = sap->get_next" << std::endl;
                //HERE
                py::tuple segTuple = sap->get_next_out_segment(true);
                // auto psegOutCap = segTuple[0].cast<py::capsule>();
                // Pstruct<Ptype::spherical_r>* psegOut = psegInCap;
                // # std::cout << "psegOut = segTuple" << std::endl;
                psegOut = segTuple[0].cast<py::capsule>(); // This compresses the above 2 lines into 1
                
                // If using Numpy structured arrays:
                /*
                py::tuple segTuple = sap->get_next_out_segment();
                auto psegOutNumpyArray = segTuple[0].cast<py::array_t<Pstruct<Ptype::spherical_r>,0>>();
                const auto psegOutNumpyArray_info = psegOutNumpyArray.request();
                const auto psegOut = static_cast<Pstruct<Ptype::spherical_r>*>(psegOutNumpyArray_info.ptr);
                */
                // # std::cout << "Reset the counter" << std::endl;                
                ipOut = 0; // Reset the counter for the new segment
              }
            
            psegOut[ipOut] = psegIn[ipIn]; // Copy all of this particle's data from the
                                           // input slot to the output slot

            // # std::cout << "Accelerate the charged" << std::endl;
            // Accelerate the charged particle with the field
            if (EsegPtr != nullptr)
              {
                auto EsegProxy = EsegPtr->mutable_unchecked<2>();
                psegOut[ipOut].ur_ = psegIn[ipIn].ur_ + qmdt*EsegProxy(ipIn, 0);
                //psegOut[ipOut].uy_ = psegIn[ipIn].uy_ + qmdt*EsegProxy(ipIn, 1);
                //psegOut[ipOut].uz_ = psegIn[ipIn].uz_ + qmdt*EsegProxy(ipIn, 2);
              }

            // Move the charged particle
            // 1. Save the initial position
            psegOut[ipOut].r0_ = psegOut[ipOut].r_;
            //psegOut[ipOut].y0_ = psegOut[ipOut].y_;
            //psegOut[ipOut].z0_ = psegOut[ipOut].z_;
            // 2. Move the particle. NON-RELATIVISTIC: v and u are the same
            psegOut[ipOut].r_ += psegOut[ipOut].ur_*dt;
            //psegOut[ipOut].y_ += psegOut[ipOut].uy_*dt;
            //psegOut[ipOut].z_ += psegOut[ipOut].uz_*dt;

            auto pCellIndex = psegOut[ipOut].cell_index_;

            //std::cout << "advance_charged_species_in_E_field. step: " << step << " ipOut: " << ipOut << " x_: " << psegOut[ipOut].r_ << " y_: " << psegOut[ipOut].y_ << std::endl;
            
            // Loop until the particle is in the current cell
            auto mLastFacet = NO_FACET;
            size_t facetCrossCount = 0;
            auto tStart = time - dt;
            // # std::cout << "dtRemaining" << std::endl;
            auto dtRemaining = dt;

            // Create a view to this cell in the dolfin::Mesh object so that we can
            // call the entities() function to get indices of vertices, facets, etc.
            // dolfin::Cell* pcellPtr = new dolfin::Cell(mesh, pCellIndex);
            pcellPtr = new dolfin::Cell(mesh, pCellIndex);
            auto vertices = pcellPtr->entities(0);

            // std::cout << "vertices " << vertices[0] << ", " << vertices[1] << std::endl;
            
            particlePosition[0] = psegOut[ipOut].r_;
            //particlePosition[1] = psegOut[ipOut].y_;
            //particlePosition[2] = psegOut[ipOut].z_;
            // std::cout << "Calling is_inside_vertices" << std::endl;

            // # std::cout << "is_inside_vertices" << std::endl;
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
                    std::cout << "DnT:particle.h:advance_charged_species_in_E_field: !!! MAX_FACET_CROSS_COUNT exceeded!!!" << std::endl;
                    exit(EXIT_FAILURE);
                  }
                
                // Compute dx[], the move vector that starts in the current cell
                // Current particle position
                pCoord2[0] = psegOut[ipOut].r_;
                // pCoord2[1] = psegOut[ipOut].y_;
                // pCoord2[2] = psegOut[ipOut].z_;
                
                // Starting particle position in this cell
                // pCoord2[3] = psegOut[ipOut].r0_;
                // pCoord2[4] = psegOut[ipOut].y0_;
                // pCoord2[5] = psegOut[ipOut].z0_;
                pCoord2[pDim] = psegOut[ipOut].r0_;
                // pCoord2[pDim+1] = psegOut[ipOut].y0_;
                // pCoord2[pDim+2] = psegOut[ipOut].z0_;

                // The move vector in this cell
                //dx[0] = pCoord2[0] - pCoord2[3];
                //dx[1] = pCoord2[1] - pCoord2[4];
                //dx[2] = pCoord2[2] - pCoord2[5];
                dx[0] = pCoord2[0] - pCoord2[pDim];
                // dx[1] = pCoord2[1] - pCoord2[pDim+1];
                // dx[2] = pCoord2[2] - pCoord2[pDim+2];

                //std::cout << "advance_charged_species_in_E_field. pCoord2 x_: " << pCoord2[pDim] << " " << pCoord2[pDim+1] << " dx " << dx[0] << " " << dx[1] << std::endl;                
                // bool returnStdArray = true;
                bool returnStdArray = false;
                // # std::cout << "facetTupl" << std::endl;                
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
                    psegOut[ipOut].r0_ = psegOut[ipOut].r0_ + dxFraction*dx[0];
                    // psegOut[ipOut].y0_ = psegOut[ipOut].y0_ + dxFraction*dx[1];
                    // psegOut[ipOut].z0_ = psegOut[ipOut].z0_ + dxFraction*dx[2];
                    
                    // Look up the mesh-level index of this facet...
                    auto mFacetArray = pcellPtr->entities(tDim-1); // The mesh-level indices of the facets of this cell.
                    auto mFacet = mFacetArray[cFacet]; // The mesh-level index of the facet just crossed.

                    // mFacet should never be the same as the last facet crossed: check this
                    if (mFacet == (size_t)mLastFacet) // If the particle has crossed the same facet twice in succession, there's an error:
                      {
                        std::string errorMsg = "In advance_charged_species_in_E_field(): The mesh index of the facet crossed is "  + std::to_string(mFacet) + ", the same as the last facet crossed. This should not happen since the particle cannot cross the same facet twice in one move!";
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
                    // # std::cout << "facValue" << std::endl;
                    if (facValue != 0)
                      {
                        // Call the function associated with this value.
                        if ((psegOut[ipOut].bitflags_ & Pstruct<Ptype::spherical_r>::TRAJECTORY_FLAG) != 0b0) // If this is a trajectory particle.
                          {
                            // Get the storage index that currently identifies this
                            // particle in the trajectory list of particles.
                            auto fullIndex = sap->get_full_index(ipIn, "in");
                            //beg
                            // Get the current solved E at the particle
                            auto EsegProxy = EsegPtr->mutable_unchecked<2>();
                            auto E_arrProxy = E1.mutable_unchecked<2>();
                            E_arrProxy(0, 0) = EsegProxy(ipIn, 0);
                            // E_arrProxy(0, 1) = EsegProxy(ipIn, 1);
                            // E_arrProxy(0, 2) = EsegProxy(ipIn, 2);
                            // External E
                            auto EextProxy = Eext.mutable_unchecked<2>();
                            auto Eext_arrProxy = Eext1.mutable_unchecked<2>();
                            Eext_arrProxy(0, 0) = EextProxy(ipIn, 0);
                            // Eext_arrProxy(0, 1) = EextProxy(ipIn, 1);
                            // Eext_arrProxy(0, 2) = EextProxy(ipIn, 2);
                            //end
                            // Record the data for the trajectory
                            particle_P.attr("record_trajectory_datum")(std::string(species_name), ipOut, fullIndex, step, tStart, "E_field_comps"_a=E1, "external_E_field_comps"_a=Eext1, "facet_crossing"_a=true); // Does this work?
                          }
                        //std::cout << "Hello from advance_charged_species_in_E_field_spherical_r@6" << std::endl;
                        // # std::cout << "advance_charged_species_in_E_field_spherical_r" << " calling a bcFunction" << std::endl;
                        auto bcFunction = bcFunctionDict[std::make_pair(facValue,species_name)];
                        (userPBFs.*bcFunction)(psegOut[ipOut], species_name, mFacet, dx, dxFraction, facetNormalVector);
                        // # std::cout << "advance_charged_species_in_E_field_spherical_r" << " returned from bcFunction" << std::endl;                        
                        
                      }

                    // Look up the cell index of the new cell.
                    // # std::cout << "pCellIndexNew" << std::endl;
                    auto pCellIndexNew = meshEntityArrays->get_cell_neighbors(pCellIndex)[cFacet];

                    /* If the particle has left the mesh, and has been deleted, end
                       the search.  If it hasn't been deleted (e.g., reflected),
                       continue tracking it.
                    */
                    if (pCellIndexNew == NO_CELL)
                      {
                        if ((psegOut[ipOut].bitflags_ & Pstruct<Ptype::spherical_r>::DELETE_FLAG) == 0b1)
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
                    std::string errorMsg = "In advance_charged_species_in_E_field_spherical_r(): The cell index of the facet crossed is NO_FACET (" + std::to_string(cFacet) + "). This should not happen since the particle has left its initial cell!";
                    std::cout << errorMsg << std::endl;                    
                    exit(EXIT_FAILURE);
                  } // END: if (cFacet != NO_FACET)

                // The following is needed in C++ because of is_inside_vertices() usage.
                // Prepare for the next 'while' iteration: need the new cell index and the cell vertices:
                delete pcellPtr; // To avoid a leak
                pcellPtr = new dolfin::Cell(mesh, pCellIndex);
                vertices = pcellPtr->entities(0);
                
                particlePosition[0] = psegOut[ipOut].r_;
                //particlePosition[1] = psegOut[ipOut].y_;
                //particlePosition[2] = psegOut[ipOut].z_;
              } // END: while (! is_inside_vertices(mesh, vertices, particlePosition))
            delete pcellPtr; // Not needed past this point.

            // Record the number of facet-crossings
            psegOut[ipOut].crossings_ = facetCrossCount;

            if ((psegOut[ipOut].bitflags_ & Pstruct<Ptype::spherical_r>::DELETE_FLAG) == 0b0) // If this particle has not been deleted...
              {
                /*
                  If particle indices in the "out" SA have changed (i.e., indexChange
                  is 'true') due to deletions, then update this particle's index where
                  needed.  E.g., if this is a trajectory particle, update its SA
                  index in the list of trajectory particles.
                */
                // # std::cout << "!traj_T.is_none" << std::endl;                
                if (!traj_T.is_none())
                  {
                    if ((psegOut[ipOut].bitflags_ & Pstruct<Ptype::spherical_r>::TRAJECTORY_FLAG) != 0b0) // If this is a trajectory particle.
                      {
                        if (indexChange == true)
                          {
                            particle_P.attr("update_trajectory_particleId")(species_name, ipIn, ipOut);
                          }
                      }
                  }
                // Advance the "out" array counter for the next particle
                // # std::cout << "ipOut += 1" << std::endl;
                ipOut += 1;
              }
            else // This particle has been deleted
              {
                indexChange = true; // This indicates that particle SA indices are changed past this point due to deletions.
                // If this was a trajectory particle, remove it's index from the trajectory-particle list.
                if (!traj_T.is_none())
                  {
                    if ((psegIn[ipIn].bitflags_ & Pstruct<Ptype::spherical_r>::TRAJECTORY_FLAG) != 0b0)
                      {
                        // # std::cout << "remove_trajectory_particleId" << std::endl;                        
                        particle_P.attr("remove_trajectory_particleId")(species_name, ipIn, ipOut, step, time, dt);
                      }
                  }
              }

            /*
              Check if we've reached the end of this segment.  If so, we need to start
              writing on a new segment.  If there are no more segments, allocate a new one.
            */
            // # std::cout << "ipOut == segmentLength" << std::endl;
            if (ipOut == segmentLength)
              {
                particleCount += segmentLength;
              }

          } // loop over particles in the "in" segment
        
        // Done with this "in" segment. Get the next one, if it exists.
        // # std::cout << "segTuple = sap->" << std::endl;
        py::tuple segTuple = sap->get_next_segment("in", true);
        auto npSeg = segTuple[0].cast<py::ssize_t>();
        if (npSeg >  0)
          {
            psegIn = segTuple[1].cast<py::capsule>(); // This is a compressed expression: there's a cast performed to make psegIn, which has it's type from above (Pstruct<Ptype::spherical_r>* psegIn)
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
        auto psegInNumpyArray = segTuple[1].cast<py::array_t<Pstruct<Ptype::spherical_r>,0>>();
        const auto psegInNumpyArray_info = psegInNumpyArray.request();
        psegIn = static_cast<Pstruct<Ptype::spherical_r>*>(psegInNumpyArray_info.ptr);
        */
        // # std::cout << "keep looping" << std::endl;
      } // while there are more "in" segments, keep looping

    // Set values that will be used to keep track of the particle arrays
    if (ipOut != segmentLength) // This catches the case where we exit the loop when ipOut = segmentLength and there are no more "in" segments.  Otherwise, we would add segmentLength to particleCount twice.
      {
        particleCount += ipOut;
      }
    //    std::cout << "particleCount: " << particleCount << std::endl;
    sap->set_number_of_items("out", particleCount);
        
  }
  // ENDDEF: void advance_charged_species_in_E_field(py::object particle_P, py::str species_name, py::object ctrl)
