  //! Advance a cartesian_xyz charged-particle species one timestep in uniform fields.
  /*!          
    Apply the electric field ctrl.E0 to particles of the given species for a
    time interval ctrl.dt.  Compute the resulting change in particle
    velocities and positions.

    NB: There is no implementation of boundary-conditions in this function, so
    particles do not get deleted.

    \param particle_P is a Particle_C object.
    \param species_name is the py::str name of the charged-particle species.
    \param ctrl is a DTcontrol_C object.

    \var SegmentedArrayPair<Ptype::cartesian_xyz>* sap
    This is a pointer-to-SegmentedArrayPair containing the particle data for one neutral species.

  */
  //  template <Ptype PT>
  //template <>
  void advance_charged_species_in_uniform_fields_cartesian_xyz(py::object particle_P, py::str species_name, py::object ctrl)
  {

    //    std::cout << "Hello from advance_charged_species_in_uniform_fields@0" << std::endl;

    // Get attributes needed from Particle_C arg
    auto qom = particle_P.attr("qom")[species_name].cast<double>();
    // auto pDim = particle_P.attr("particle_dimension").cast<int>();

    // Get attributes needed from ctrl arg
    auto dt = ctrl.attr("dt").cast<double>(); // This does a COPY.
    auto E0x = ctrl.attr("E0").attr("x").cast<double>();
    auto E0y = ctrl.attr("E0").attr("y").cast<double>();
    auto E0z = ctrl.attr("E0").attr("z").cast<double>();

    auto qmdt = qom*dt;

    //    std::cout << "pDim " << pDim << std::endl;

    //    std::cout << "Hello from advance_charged_species_in_uniform_fields@1" << std::endl;

    // Get the data for the species to be advanced
    // auto sap_map = particle_P.attr("sap_dict").cast<std::map<std::string, SegmentedArrayPair<Ptype::cartesian_xyz> *>>();
    // auto sap = sap_map[std::string(species_name)];
    auto sap = particle_P.attr("sap_dict")[species_name].cast<SegmentedArrayPair<Ptype::cartesian_xyz> *>();
    
    const py::ssize_t segmentLength = sap->get_segment_length();
    
    // Start a loop over the sap segments. segTuple contains (npSeg, psegIn, psegOut)
    py::tuple segTuple = sap->init_inout_loop(true); // 'true' returns the data pointers for
                                                     // the first segment pair, instead of the
                                                     // Numpy arrays.
    auto npSeg = segTuple[0].cast<py::ssize_t>();

    // Get the particle array pointers from the returned py::capsules.
    // psegIn:
    auto psegInCap = segTuple[1].cast<py::capsule>();
    // std::cout << "The name of the capsule is " << psegInCap.name() << std::endl;
    Pstruct<Ptype::cartesian_xyz>* psegIn = psegInCap; // Declares the type of psegIn
    // std::cout << "psegIn[0].x_ = " << psegIn[0].x_ << std::endl;
    // psegOut:
    auto psegOutCap = segTuple[2].cast<py::capsule>();
    // std::cout << "The name of the capsule is " << psegOutCap.name() << std::endl;
    Pstruct<Ptype::cartesian_xyz>* psegOut = psegOutCap; // Declares the type of psegOut

    // An alternative way to do this is to use the returned Numpy structured array
    // objects: call init_inout_loop() without the "true" argument above. Then use the
    // Python buffer protocol to get at the data.
    /*
    // psegIn:
    auto psegInNumpyArray = segTuple[1].cast<py::array_t<Pstruct<Ptype::cartesian_xyz>,0>>();
    const auto psegInNumpyArray_info = psegInNumpyArray.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
    auto psegIn = static_cast<Pstruct<Ptype::cartesian_xyz>*>(psegInNumpyArray_info.ptr); // Pointer to a particle struct in psegIn
    // psegOut:
    auto psegOutNumpyArray = segTuple[2].cast<py::array_t<Pstruct<Ptype::cartesian_xyz>,0>>();
    const auto psegOutNumpyArray_info = psegOutNumpyArray.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
    auto psegOut = static_cast<Pstruct<Ptype::cartesian_xyz>*>(psegOutNumpyArray_info.ptr); // Pointer to a particle struct in psegOut
    */
    
    py::ssize_t particleCount = 0; // Counts the number of still-active particles for this
                                   // species.
    
    // ipOut counts particles being written to the current "out" segment.
    py::ssize_t ipOut = 0;
    // bool indexChange = false; // Flag to indicate if the particle SA indices have
                              // changed. This affects, e.g., trajectories, which use SA
                              // indices to identify particles chosen for trajectory
                              // plots.

    //    std::cout << "Hello from advance_charged_species_in_uniform_fields@2" << std::endl;
    
    while (psegIn != nullptr) // Keep looping until we run out of "in" segments.
      {
        for (auto ipIn = 0; ipIn < npSeg; ipIn++)
          {
            /* No particle deletion in this algorithm
            // Skip deleted particles
            if ((psegIn[ipIn].bitflags_ & Pstruct<Ptype::cartesian_xyz>::DELETE_FLAG) != 0b0)
              {
                indexChange = true; // Particle SA indices are stale past this point
                                    // due to deletions.
                continue; // skip to the next particle in the "in" segment.
              }
            */

            // If this "out" segment is full, get the next "out" segment, since
            // there are still some "in" particles to advance.  If there are no more
            // "out" segments, allocate a new one.
            if (ipOut == segmentLength)
              {
                // Using pointers to the plain arrays:
                py::tuple segTuple = sap->get_next_out_segment(true);                
                // auto psegOutCap = segTuple[0].cast<py::capsule>();
                // Pstruct<Ptype::cartesian_xyz>* psegOut = psegInCap;
                psegOut = segTuple[0].cast<py::capsule>(); // This compresses the above 2 lines into 1
                
                // Using Numpy structured arrays instead:
                /*
                py::tuple segTuple = sap->get_next_out_segment();
                auto psegOutNumpyArray = segTuple[0].cast<py::array_t<Pstruct<Ptype::cartesian_xyz>,0>>();
                const auto psegOutNumpyArray_info = psegOutNumpyArray.request();
                const auto psegOut = static_cast<Pstruct<Ptype::cartesian_xyz>*>(psegOutNumpyArray_info.ptr);
                */
                
                ipOut = 0; // Reset the counter for the new segment
              }
            
            psegOut[ipOut] = psegIn[ipIn]; // Copy all of this particle's data from the
                                           // input slot to the output slot

            // Accelerate the particle
            psegOut[ipOut].ux_ = psegIn[ipIn].ux_ + qmdt*E0x;
            psegOut[ipOut].uy_ = psegIn[ipIn].uy_ + qmdt*E0y;
            psegOut[ipOut].uz_ = psegIn[ipIn].uz_ + qmdt*E0z;

            // Move the particle
            // 1. Save the initial position (No mesh, so not needed?)
            psegOut[ipOut].x0_ = psegOut[ipOut].x_;
            psegOut[ipOut].y0_ = psegOut[ipOut].y_;
            psegOut[ipOut].z0_ = psegOut[ipOut].z_;
            // 2. Move the particle. NON-RELATIVISTIC: v and u are the same
            psegOut[ipOut].x_ += psegOut[ipOut].ux_*dt;
            psegOut[ipOut].y_ += psegOut[ipOut].uy_*dt;
            psegOut[ipOut].z_ += psegOut[ipOut].uz_*dt;

            ipOut += 1;
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
        py::tuple segTuple = sap->get_next_segment("in");
        npSeg = segTuple[0].cast<py::ssize_t>();
        auto psegInNumpyArray = segTuple[1].cast<py::array_t<Pstruct<Ptype::cartesian_xyz>,0>>();
        const auto psegInNumpyArray_info = psegInNumpyArray.request();
        psegIn = static_cast<Pstruct<Ptype::cartesian_xyz>*>(psegInNumpyArray_info.ptr);
        */
      } // while there are more "in" segments, keep looping

    // Set values that will be used to keep track of the particle arrays
    if (ipOut != segmentLength) // This catches the case where we exit the loop when ipOut = segmentLength and there are no more "in" segments.  Otherwise, we would add segmentLength to particleCount twice.
      {
        particleCount += ipOut;
      }
    sap->set_number_of_items("out", particleCount);
  }
  // ENDDEF: void advance_charged_species_in_uniform_fields(py::object particle_P, py::str species_name, py::object ctrl)
