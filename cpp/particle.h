/*! \file particle.h

  \brief This file has the source code for a C++ implementation of particle advance.

  \namespace dnt

  template <Ptype PT>
    void advance_neutral_species(py::object particle_P, py::str species_name, py::object ctrl)

  template <Ptype PT>
    void advance_charged_species_in_uniform_fields(py::object particle_P, py::str species_name, py::object ctrl)

  \sa particle_solib.cpp, MeshEntityArrays.h, SegmentedArrayPair.h, ParticleMeshBoundaryConditions.h, UserParticleBoundaryFunctions.h, user_particle_boundary_functions_solib.cpp
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
#include "SegmentedArrayPair.h"
#include "MeshEntityArrays.h"
#include "ParticleMeshBoundaryConditions.h"
#include "UserParticleBoundaryFunctions.h"

#include "dolfin_functions.h"

namespace py = pybind11;

// Static scratch space
static double particlePosition[3];
// A scratch array to hold: x,y,z, x0,y0,z0 (or a subset)
static double pCoord2[6];
// A scratch array to hold: dx,dy,dz (or a subset)
static double dx[3];

static dolfin::Cell* pcellPtr; // File scoped needed? Probably not.

// From https://en.cppreference.com/w/cpp/language/function_template:
// A function template by itself is not a type, or a function, or any other entity. No
// code is generated from a source file that contains only template definitions. In order
// for any code to appear, a template must be instantiated: the template arguments must be
// determined so that the compiler can generate an actual function (or class, from a class
// template).

namespace dnt
{

    //! Create objects needed to advance a particle species.
    /*!
      The Numpy arrays passed in have been allocated in Particle_Module.py.
      
      \return void

    */
  void initialize_particle_integration(py::array_t<double> negE_in, py::array_t<double> Eext_in, py::array_t<double> zeroE_in)
    {
      std::cout << "Hello from initialize_particle_integration@1" << std::endl;
      // Initialize the static interpolation arrays

      py::buffer_info negEinfo = negE_in.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
      negE_array = (double*) negEinfo.ptr;
      E_array_len = negEinfo.shape[0];
      E_array_ncomps = negEinfo.shape[1];

      py::buffer_info EextInfo = Eext_in.request();
      Eext_array = (double*) EextInfo.ptr;
      // Check that the arrays have the same size
      if ((EextInfo.shape[0] != E_array_len) || (EextInfo.shape[1] != E_array_ncomps))
        throw std::runtime_error("particle.h:initialize_particle_integration@1: Array dimensions are not the same");

      py::buffer_info zeroEinfo = zeroE_in.request();
      zeroE_array = (double*) zeroEinfo.ptr;
      // Check that the arrays have the same size
      if ((zeroEinfo.shape[0] != E_array_len) || (zeroEinfo.shape[1] != E_array_ncomps))
        throw std::runtime_error("particle.h:initialize_particle_integration@2: Array dimensions are not the same");
      */
    }
  // ENDDEF: void initialize_particle_integration
  
  
  //! Advance a charged-particle species by one time-increment on a mesh.
  /*!          

    Apply interpolated electric force to charged particles.  Compute the change in
    velocity and position in time dt. Use an explicit method to integrate the orbit.

    If a particle leaves its initial cell, the cell that the particle moves to is
    calculated by finding what facet the particle crosses, and looking up what the
    neighbor cell is. This is repeated until the cell containing the final position
    is found.

    \param particle_P is a Particle_C Python object.
    \param species_name is the py::str name of the species to be advanced.
    \param ctrl is a DTcontrol_C Python object.

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
  // Make a special for cartesian_xy, triangle mesh:

//  template <>
// Make N_CELL_FACETS the only template parameter:
  template <size_t N_CELL_FACETS>
    void advance_charged_species_in_E_field_cartesian_xy(py::object &particle_P, py::str species_name, py::object &ctrl, dolfin::Function *neg_E_field, dolfin::Function *external_E_field, bool accel_only)
  {
    using namespace pybind11::literals;

    //std::cout << "Hello from advance_charged_species_in_E_field_cartesian_xy@1" << std::endl;
    //std::cout << "This is species " << std::string(species_name) << std::endl;

    // Maybe all of the stuff that doesn't change from call to call should be moved to "initialize", and used to initialize static variables?
    // Marked these with (m) below
    // For parallel runs, need to call "intialize" again after a rebalance.
    
    // Get attributes from the Particle_C argument
    
    auto pmesh_M = particle_P.attr("pmesh_M"); // (m)
    auto pDim = particle_P.attr("particle_dimension").cast<int>(); // (m)
    auto sap = particle_P.attr("sap_dict")[species_name].cast<SegmentedArrayPair<Ptype::PARTICLE_TYPE> *>();
    // Note that sap is a pointer.
    // Could also get sap by first casting to an std::map:
    //   auto sap_map = particle_P.attr("sap_dict").cast<std::map<std::string, SegmentedArrayPair<PT> *>>();
    //   auto sap = sap_map[std::string(species_name)];
    //   auto sapMap = particle_P.attr("sap_dict").cast<std::map<std::string, SegmentedArrayPair<Ptype::cartesian_xyz> *>>();
    auto qom = particle_P.attr("qom")[species_name].cast<double>();
    // Get call-back functions for boundary-crossing particles
    // Move this to the top of the function
    // This was when we were using the Python callbacks
    // auto bcFunctionDict = particle_P.attr("pmesh_bcs").attr("bc_function_dict").cast<std::map<int, py::dict>>();
    auto pMBCs = particle_P.attr("pmesh_bcs").cast<ParticleMeshBoundaryConditions<Ptype::PARTICLE_TYPE>>(); // (m)
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
        //std::cout << "advance_charged_species_in_E_field_cartesian_xy: applySolvedElectricFieldMap[species_name] = " << applySolvedElectricFieldMap[species_name] << std::endl;
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
        //std::cout << "advance_charged_species_in_E_field_cartesian_xy: applyRandomExternalElectricFieldMap[species_name] = " << applyRandomExternalElectricFieldMap[species_name] << std::endl;
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
    Pstruct<Ptype::PARTICLE_TYPE>* psegIn = psegInCap; // Declares the type of psegIn
    // std::cout << "psegIn[0].x_ = " << psegIn[0].x_ << std::endl;
    // psegOut:
    auto psegOutCap = segTuple[2].cast<py::capsule>();
    // std::cout << "The name of the capsule is " << psegOutCap.name() << std::endl;
    Pstruct<Ptype::PARTICLE_TYPE>* psegOut = psegOutCap; // Declares the type of psegOut

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

    //std::cout << "Hello from advance_charged_species_in_E_field_cartesian_xy@2" << std::endl;
    
    py::array_t<double>* EsegPtr = nullptr;
    while (psegIn != nullptr) // Keep looping until we run out of "in" segments.
      {

        //std::cout << "Hello from advance_charged_species_in_E_field_cartesian_xy@3" << std::endl;
        // Compute the electric field at each particle

        // Eseg = None # This is for the case where no fields are to be applied to the particles,
        //        if (!neg_E_field.is_none())
        if (neg_E_field != nullptr)
          {
            //std::cout << "Hello from advance_charged_species_in_E_field_cartesian_xy@3.1" << " applySolvedElectricFieldFlag = " << applySolvedElectricFieldFlag << std::endl;
            if ((applySolvedElectricFieldFlag == true) || (applySolvedElectricFieldMap[species_name] == true))
              //if (applySolvedElectricFieldFlag == true)
              {
                //std::cout << "Hello from advance_charged_species_in_E_field_cartesian_xy@3.2" << " npSeg = " << npSeg << std::endl;
                interpolate_field_to_points<Ptype::PARTICLE_TYPE>(neg_E_field, psegIn, npSeg, negE); // From dolfin_functions.h
                //std::cout << "Hello from advance_charged_species_in_E_field_cartesian_xy@4" << std::endl;
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
            //std::cout << "Hello from advance_charged_species_in_E_field_cartesian_xy@5" << std::endl;
            EsegPtr = &zeroE;
          }
        // Add an external E field, if present.
        if (external_E_field != nullptr)
          {
            if ((applyRandomExternalElectricFieldFlag == true) || (applyRandomExternalElectricFieldMap[species_name] == true))
              {
                // Interpolate the external electric field to the particle positions
                std::cout << "advance_charged_species_in_E_field_cartesian_xy: This function needs to be written" << std::endl;                
                // interpolate_random_field_to_points<Ptype::PARTICLE_TYPE>(external_E_field, psegIn, npSeg, Eext);
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
            if (ipOut == segmentLength)
              {
                // Using pointers to the plain arrays:
                py::tuple segTuple = sap->get_next_out_segment(true);                
                // auto psegOutCap = segTuple[0].cast<py::capsule>();
                // Pstruct<Ptype::PARTICLE_TYPE>* psegOut = psegInCap;
                psegOut = segTuple[0].cast<py::capsule>(); // This compresses the above 2 lines into 1
                
                // If using Numpy structured arrays:
                /*
                py::tuple segTuple = sap->get_next_out_segment();
                auto psegOutNumpyArray = segTuple[0].cast<py::array_t<Pstruct<Ptype::PARTICLE_TYPE>,0>>();
                const auto psegOutNumpyArray_info = psegOutNumpyArray.request();
                const auto psegOut = static_cast<Pstruct<Ptype::PARTICLE_TYPE>*>(psegOutNumpyArray_info.ptr);
                */
                
                ipOut = 0; // Reset the counter for the new segment
              }
            
            psegOut[ipOut] = psegIn[ipIn]; // Copy all of this particle's data from the
                                           // input slot to the output slot

            // Accelerate the charged particle with the field
            if (EsegPtr != nullptr)
              {
                auto EsegProxy = EsegPtr->mutable_unchecked<2>();
                psegOut[ipOut].ux_ = psegIn[ipIn].ux_ + qmdt*EsegProxy(ipIn, 0);
                psegOut[ipOut].uy_ = psegIn[ipIn].uy_ + qmdt*EsegProxy(ipIn, 1);
                //psegOut[ipOut].uz_ = psegIn[ipIn].uz_ + qmdt*EsegProxy(ipIn, 2);
              }

            // Move the charged particle
            // 1. Save the initial position
            psegOut[ipOut].x0_ = psegOut[ipOut].x_;
            psegOut[ipOut].y0_ = psegOut[ipOut].y_;
            //psegOut[ipOut].z0_ = psegOut[ipOut].z_;
            // 2. Move the particle. NON-RELATIVISTIC: v and u are the same
            psegOut[ipOut].x_ += psegOut[ipOut].ux_*dt;
            psegOut[ipOut].y_ += psegOut[ipOut].uy_*dt;
            //psegOut[ipOut].z_ += psegOut[ipOut].uz_*dt;

            auto pCellIndex = psegOut[ipOut].cell_index_;

            //std::cout << "advance_charged_species_in_E_field. step: " << step << " ipOut: " << ipOut << " x_: " << psegOut[ipOut].x_ << " y_: " << psegOut[ipOut].y_ << std::endl;
            
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
            //particlePosition[2] = psegOut[ipOut].z_;
            // std::cout << "Calling is_inside_vertices" << std::endl;
            
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
                pCoord2[0] = psegOut[ipOut].x_;
                pCoord2[1] = psegOut[ipOut].y_;
                // pCoord2[2] = psegOut[ipOut].z_;
                
                // Starting particle position in this cell
                // pCoord2[3] = psegOut[ipOut].x0_;
                // pCoord2[4] = psegOut[ipOut].y0_;
                // pCoord2[5] = psegOut[ipOut].z0_;
                pCoord2[pDim] = psegOut[ipOut].x0_;
                pCoord2[pDim+1] = psegOut[ipOut].y0_;
                // pCoord2[pDim+2] = psegOut[ipOut].z0_;

                // The move vector in this cell
                //dx[0] = pCoord2[0] - pCoord2[3];
                //dx[1] = pCoord2[1] - pCoord2[4];
                //dx[2] = pCoord2[2] - pCoord2[5];
                dx[0] = pCoord2[0] - pCoord2[pDim];
                dx[1] = pCoord2[1] - pCoord2[pDim+1];
                // dx[2] = pCoord2[2] - pCoord2[pDim+2];

                //std::cout << "advance_charged_species_in_E_field. pCoord2 x_: " << pCoord2[pDim] << " " << pCoord2[pDim+1] << " dx " << dx[0] << " " << dx[1] << std::endl;                
                // bool returnDataPtr = true;
                bool returnDataPtr = false;
                py::tuple facetTupl = find_facet(pmesh_M, &pCoord2[pDim], dx, pCellIndex, returnDataPtr);
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
                    //psegOut[ipOut].z0_ = psegOut[ipOut].z0_ + dxFraction*dx[2];
                    
                    // Look up the mesh-level index of this facet...
                    auto mFacetArray = pcellPtr->entities(tDim-1); // The mesh-level indices of the facets of this cell.
                    auto mFacet = mFacetArray[cFacet]; // The mesh-level index of the facet just crossed.

                    // mFacet should never be the same as the last facet crossed: check this
                    if (mFacet == (size_t)mLastFacet) // If the particle has crossed the same facet twice in succession, there's an error:
                      {
                        std::string errorMsg = "In advance_charged_species_in_E_field(): The mesh index of the facet crossed is "  + std::to_string(mFacet) + ", the same as the last facet crossed. This should not happen since the particle cannot cross the same facet twice in one move!";
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
                        if ((psegOut[ipOut].bitflags_ & Pstruct<Ptype::PARTICLE_TYPE>::TRAJECTORY_FLAG) != 0b0) // If this is a trajectory particle.
                          {
                            // Get the storage index that currently identifies this
                            // particle in the trajectory list of particles.
                            auto fullIndex = sap->get_full_index(ipIn, "in");
                            //beg
                            // Get the current solved E at the particle
                            auto EsegProxy = EsegPtr->mutable_unchecked<2>();
                            auto E_arrProxy = E1.mutable_unchecked<2>();
                            E_arrProxy(0, 0) = EsegProxy(ipIn, 0);
                            E_arrProxy(0, 1) = EsegProxy(ipIn, 1);
                            // E_arrProxy(0, 2) = EsegProxy(ipIn, 2);
                            // External E
                            auto EextProxy = Eext.mutable_unchecked<2>();
                            auto Eext_arrProxy = Eext1.mutable_unchecked<2>();
                            Eext_arrProxy(0, 0) = EextProxy(ipIn, 0);
                            Eext_arrProxy(0, 1) = EextProxy(ipIn, 1);
                            // Eext_arrProxy(0, 2) = EextProxy(ipIn, 2);
                            //end
                            // Record the data for the trajectory
                            particle_P.attr("record_trajectory_datum")(std::string(species_name), ipOut, fullIndex, step, tStart, "E_field_comps"_a=E1, "external_E_field_comps"_a=Eext1, "facet_crossing"_a=true); // Does this work?
                          }
                        //std::cout << "Hello from advance_charged_species_in_E_field_cartesian_xy@6" << std::endl;
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
                        if ((psegOut[ipOut].bitflags_ & Pstruct<Ptype::PARTICLE_TYPE>::DELETE_FLAG) == 0b1)
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
                    exit(EXIT_FAILURE);
                  } // END: if (cFacet != NO_FACET)

                // The following is needed in C++ because of is_inside_vertices() usage.
                // Prepare for the next 'while' iteration: need the new cell index and the cell vertices:
                delete pcellPtr; // To avoid a leak
                pcellPtr = new dolfin::Cell(mesh, pCellIndex);
                vertices = pcellPtr->entities(0);
                
                particlePosition[0] = psegOut[ipOut].x_;
                particlePosition[1] = psegOut[ipOut].y_;
                //particlePosition[2] = psegOut[ipOut].z_;
              } // END: while (! is_inside_vertices(mesh, vertices, particlePosition))
            delete pcellPtr; // Not needed past this point.

            // Record the number of facet-crossings
            psegOut[ipOut].crossings_ = facetCrossCount;

            if ((psegOut[ipOut].bitflags_ & Pstruct<Ptype::PARTICLE_TYPE>::DELETE_FLAG) == 0b0) // If this particle has not been deleted...
              {
                /*
                  If particle indices in the "out" SA have changed (i.e., indexChange
                  is 'true') due to deletions, then update this particle's index where
                  needed.  E.g., if this is a trajectory particle, update its SA
                  index in the list of trajectory particles.
                */
                if (!traj_T.is_none())
                  {
                    if ((psegOut[ipOut].bitflags_ & Pstruct<Ptype::PARTICLE_TYPE>::TRAJECTORY_FLAG) != 0b0) // If this is a trajectory particle.
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
                    if ((psegIn[ipIn].bitflags_ & Pstruct<Ptype::PARTICLE_TYPE>::TRAJECTORY_FLAG) != 0b0)
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
        auto npSeg = segTuple[0].cast<py::ssize_t>();
        if (npSeg >  0)
          {
            psegIn = segTuple[1].cast<py::capsule>(); // This is a compressed expression: there's a cast performed to make psegIn, which has it's type from above (Pstruct<Ptype::PARTICLE_TYPE>* psegIn)
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
        auto psegInNumpyArray = segTuple[1].cast<py::array_t<Pstruct<Ptype::PARTICLE_TYPE>,0>>();
        const auto psegInNumpyArray_info = psegInNumpyArray.request();
        psegIn = static_cast<Pstruct<Ptype::PARTICLE_TYPE>*>(psegInNumpyArray_info.ptr);
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
  // ENDDEF: void advance_charged_species_in_E_field(py::object particle_P, py::str species_name, py::object ctrl)

  //! Advance a neutral particle species by one time-increment on a mesh.
  /*!          

    Compute change in position in time dt. Use an explicit method to calculate the
    final position.  If a particle leaves its initial cell, the cell that the
    particle moves to is calculated by finding what facet the particle crosses, and
    looking up what the neighbor cell is. This is repeated until the cell
    containing the final position is found.

    \param particle_P is a Particle_C object.
    \param species_name is the py::str name of the species to be advanced.
    \param ctrl is a DTcontrol_C object.

    \var SegmentedArrayPair<PT>* sap
    This is a pointer-to-SegmentedArrayPair containing the particle data for one neutral species.
    \var double dtRemaining
    The time a particle has left to move in the new cell the particle has entered.
    \var int pDim
    Number of spatial coordinates in the particle location.
    \var double tStart
    The time at which a particle starts its move in the current cell.

  */
  template <size_t N_CELL_FACETS>    
    void advance_neutral_species_cartesian_xyz(py::object particle_P, py::str species_name, py::object ctrl)
    {
    using namespace pybind11::literals;
        
    // std::cout << "Hello from advance_neutral_species@0" << std::endl;
    // std::cout << "This is species " << std::string(species_name) << std::endl;

    // Get attributes from the Particle_C argument
    auto pmesh_M = particle_P.attr("pmesh_M");
    auto pDim = particle_P.attr("particle_dimension").cast<int>();
    auto sap = particle_P.attr("sap_dict")[species_name].cast<SegmentedArrayPair<Ptype::cartesian_xyz> *>();
    // Note that sap is a pointer.
    // Could also get sap by first casting to an std::map:
    //   auto sap_map = particle_P.attr("sap_dict").cast<std::map<std::string, SegmentedArrayPair<Ptype::cartesian_xyz> *>>();
    //   auto sap = sap_map[std::string(species_name)];
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
    Pstruct<Ptype::cartesian_xyz>* psegIn = psegInCap; // Declares the type of psegIn
    //    std::cout << "psegIn[0].x_ = " << psegIn[0].x_ << std::endl;
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
    auto psegOut = static_cast<Pstruct<Ptype::cartesian_xyz>*>(psegOutNumpyArray_info.ptr); // Pointer to the array of particle structs.
    */
    
    py::ssize_t particleCount = 0; // Counts the number of still-active particles for this
                                   // species.
    
    // ipOut counts particles being written to the current "out" segment.
    py::ssize_t ipOut = 0;
    bool indexChange = false; // Flag to indicate if the particle SA indices have
                              // changed. This affects, e.g., trajectories, which use SA
                              // indices to identify particles chosen for trajectory
                              // plots.

    //    std::cout << "Hello from advance_neutral_species@2" << std::endl;
    
    while (psegIn != nullptr) // Keep looping until we run out of "in" segments.
      {
        for (auto ipIn = 0; ipIn < npSeg; ipIn++)
          {
            // Skip deleted particles
            if ((psegIn[ipIn].bitflags_ & Pstruct<Ptype::cartesian_xyz>::DELETE_FLAG) != 0b0)
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
                // Pstruct<Ptype::cartesian_xyz>* psegOut = psegInCap;
                psegOut = segTuple[0].cast<py::capsule>(); // This compresses the above 2 lines into 1
                
                // If using Numpy structured arrays:
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

            // Move the neutral particle
            // 1. Save the initial position
            psegOut[ipOut].x0_ = psegOut[ipOut].x_;
            psegOut[ipOut].y0_ = psegOut[ipOut].y_;
            psegOut[ipOut].z0_ = psegOut[ipOut].z_;
            // 2. Move the particle. NON-RELATIVISTIC: v and u are the same
            psegOut[ipOut].x_ += psegOut[ipOut].ux_*dt;
            psegOut[ipOut].y_ += psegOut[ipOut].uy_*dt;
            psegOut[ipOut].z_ += psegOut[ipOut].uz_*dt;

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
            particlePosition[2] = psegOut[ipOut].z_;
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
                pCoord2[2] = psegOut[ipOut].z_;
                // Starting particle position in this cell
                pCoord2[3] = psegOut[ipOut].x0_;
                pCoord2[4] = psegOut[ipOut].y0_;
                pCoord2[5] = psegOut[ipOut].z0_;
                // The move vector in this cell
                dx[0] = pCoord2[0] - pCoord2[3];
                dx[1] = pCoord2[1] - pCoord2[4];
                dx[2] = pCoord2[2] - pCoord2[5];

                // bool returnDataPtr = true;
                bool returnDataPtr = false;
                py::tuple facetTupl = find_facet(pmesh_M, &pCoord2[pDim], dx, pCellIndex, returnDataPtr);
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
                    psegOut[ipOut].z0_ = psegOut[ipOut].z0_ + dxFraction*dx[2];
                    
                    // Look up the mesh-level index of this facet...
                    auto mFacetArray = pcellPtr->entities(tDim-1); // The mesh-level indices of the facets of this cell.
                    auto mFacet = mFacetArray[cFacet]; // The mesh-level index of the facet just crossed.

                    // mFacet should never be the same as the last facet crossed: check this
                    if (mFacet == (size_t)mLastFacet) // If the particle has crossed the same facet twice in succession, there's an error:
                      {
                        std::string errorMsg = "In advance_neutral_species(): The mesh index of the facet crossed is "  + std::to_string(mFacet) + ", the same as the last facet crossed. This should not happen since the particle cannot cross the same facet twice in one move!";
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
                        if ((psegOut[ipOut].bitflags_ & Pstruct<Ptype::cartesian_xyz>::TRAJECTORY_FLAG) != 0b0) // If this is a trajectory particle.
                          {
                            // Get the storage index that currently identifies this
                            // particle in the trajectory list of particles.
                            auto fullIndex = sap->get_full_index(ipIn, "in");
                            particle_P.attr("record_trajectory_datum")(std::string(species_name), ipOut, fullIndex, step, tStart, "facet_crossing"_a=true); // Does this work?
                          }
                        // A reference to dx[] is available in the BC function class.
                        // ?? Test for is_none() first?
                        auto bcFunctionDict = particle_P.attr("pmesh_bcs").attr("bc_function_dict").cast<py::dict>();
                        // A py::dict needs a py::str key value. Convert facValue to a string first.
                        auto bcFunction = bcFunctionDict[py::str(std::to_string(facValue))].cast<py::dict>();
                        // bcFunction[py::str(species_name)](ipOut, species_name, mFacet, "dx_fraction"_a=dxFraction, "facet_normal"_a=facetNormalVector);                        
                        bcFunction[species_name](ipOut, species_name, mFacet, "dx_fraction"_a=dxFraction, "facet_normal"_a=facetNormalVector);                        
                      }

                    // Look up the cell index of the new cell.
                    auto pCellIndexNew = meshEntityArrays->get_cell_neighbors(pCellIndex)[cFacet];

                    /* If the particle has left the mesh, and has been deleted, end
                       the search.  If it hasn't been deleted (e.g., reflected),
                       continue tracking it.
                    */
                    if (pCellIndexNew == NO_CELL)
                      {
                        if ((psegOut[ipOut].bitflags_ & Pstruct<Ptype::cartesian_xyz>::DELETE_FLAG) == 0b1)
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
                    exit(EXIT_FAILURE);
                  } // END: if (cFacet != NO_FACET)

                // The following is needed in C++ because of is_inside_vertices() usage.
                // Prepare for the next 'while' iteration: need the new cell index and the cell vertices:
                delete pcellPtr; // To avoid a leak
                pcellPtr = new dolfin::Cell(mesh, pCellIndex);
                vertices = pcellPtr->entities(0);
                
                particlePosition[0] = psegOut[ipOut].x_;
                particlePosition[1] = psegOut[ipOut].y_;
                particlePosition[2] = psegOut[ipOut].z_;
              } // END: while (! is_inside_vertices(mesh, vertices, particlePosition))
            delete pcellPtr; // Not needed past this point.

            // Record the number of facet-crossings
            psegOut[ipOut].crossings_ = facetCrossCount;

            if ((psegOut[ipOut].bitflags_ & Pstruct<Ptype::cartesian_xyz>::DELETE_FLAG) == 0b0) // If this particle has not been deleted...
              {
                /*
                  If particle indices in the "out" SA have changed (i.e., indexChange
                  is 'true') due to deletions, then update this particle's index where
                  needed.  E.g., if this is a trajectory particle, update its SA
                  index in the list of trajectory particles.
                */
                if (!traj_T.is_none())
                  {
                    if ((psegOut[ipOut].bitflags_ & Pstruct<Ptype::cartesian_xyz>::TRAJECTORY_FLAG) != 0b0) // If this is a trajectory particle.
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
                    if ((psegIn[ipIn].bitflags_ & Pstruct<Ptype::cartesian_xyz>::TRAJECTORY_FLAG) != 0b0)
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
        auto npSeg = segTuple[0].cast<py::ssize_t>();
        if (npSeg >  0)
          {
            psegIn = segTuple[1].cast<py::capsule>(); // This is a compressed expression: there's a cast performed to make psegIn, which has it's type from above (Pstruct<Ptype::cartesian_xyz>* psegIn)
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
    //    std::cout << "particleCount: " << particleCount << std::endl;
    sap->set_number_of_items("out", particleCount);
  }
  // ENDDEF: void advance_neutral_species(py::object particle_P, py::str species_name, py::object ctrl)
  
  //  template <size_t N_CELL_FACETS>
  //    void advance_neutral_species<Ptype::cartesian_xy>(py::object particle_P, py::str species_name, py::object ctrl)
  //  template <Ptype cartesian_xy, size_t N_CELL_FACETS>
  //  template <>
  template <size_t N_CELL_FACETS>
    void advance_neutral_species_cartesian_xy(py::object particle_P, py::str species_name, py::object ctrl)
  {
    using namespace pybind11::literals;
    //using PT = Ptype::cartesian_xy; not a type
        
    // std::cout << "Hello from advance_neutral_species@0" << std::endl;
    // std::cout << "This is species " << std::string(species_name) << std::endl;

    // Get attributes from the Particle_C argument
    auto pmesh_M = particle_P.attr("pmesh_M");
    auto pDim = particle_P.attr("particle_dimension").cast<int>();
    auto sap = particle_P.attr("sap_dict")[species_name].cast<SegmentedArrayPair<Ptype::cartesian_xy> *>();
    // Note that sap is a pointer.
    // Could also get sap by first casting to an std::map:
    //   auto sap_map = particle_P.attr("sap_dict").cast<std::map<std::string, SegmentedArrayPair<Ptype::cartesian_xy> *>>();
    //   auto sap = sap_map[std::string(species_name)];
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

    //    std::cout << "Hello from advance_neutral_species@2" << std::endl;
    
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

                // bool returnDataPtr = true;
                bool returnDataPtr = false;
                py::tuple facetTupl = find_facet(pmesh_M, &pCoord2[pDim], dx, pCellIndex, returnDataPtr);
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
                        // A reference to dx[] is available in the BC function class.
                        // ?? Test for is_none() first?
                        auto bcFunctionDict = particle_P.attr("pmesh_bcs").attr("bc_function_dict").cast<py::dict>();
                        // A py::dict needs a py::str key value. Convert facValue to a string first.
                        auto bcFunction = bcFunctionDict[py::str(std::to_string(facValue))].cast<py::dict>();
                        // bcFunction[py::str(species_name)](ipOut, species_name, mFacet, "dx_fraction"_a=dxFraction, "facet_normal"_a=facetNormalVector);                        
                        bcFunction[species_name](ipOut, species_name, mFacet, "dx_fraction"_a=dxFraction, "facet_normal"_a=facetNormalVector);
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
        auto npSeg = segTuple[0].cast<py::ssize_t>();
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
  
  //! Advance a charged-particle species one timestep in uniform fields.
  /*!          
    Apply the electric field ctrl.E0 to particles of the given species for a
    time interval ctrl.dt.  Compute the resulting change in particle
    velocities and positions.

    NB: There is no implementation of boundary-conditions in this function, so
    particles do not get deleted.

    \param particle_P is a Particle_C object.
    \param species_name is the py::str name of the charged-particle species.
    \param ctrl is a DTcontrol_C object.

    \var SegmentedArrayPair<Ptype::PARTICLE_TYPE>* sap
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
  
} // namespace dnt


// Do these cause the compiler to make the specialized classes?
// The .def bindings in particle_solib.cpp seem to be sufficient to get compilation.
// That's probably because the definitions are included here in particle.h.

//template void dnt::advance_neutral_species<dnt::Ptype::cartesian_xyz, 2>(py::object,
//                                                                         py::str,
//                                                                         py::object);
//template void dnt::advance_neutral_species<dnt::Ptype::cartesian_xy, 3>(py::object,
//                                                                        py::str,
//                                                                        py::object);

// Write a helper function to make specialized ParticleMeshBoundaryConditions
// classes. This is in particle.h, reflecting the fact that the Python version is in
// Particle_Module.py. This also avoids putting a long piece of code into
// particle_solib.cpp, and avoids having to make a separate module for it.  Also,
// particle.h already has the pybind11 headers.

namespace dnt {

  // The anonymous namespace limits the scope of the functions in it to this file.
  namespace {

    //! Make ParticleMeshBoundaryConditions classes to contain UserParticleBoundaryFunctions for different Ptypes.
    /*!

      makeParticleMeshBoundaryConditions() is a 'helper' function. It creates the Python
      bindings for a specialized ParticleMeshBoundaryConditions class and its member
      functions based on the PT template parameter. Class instances can then be created
      and used in Python.

      \param PT is a template parameter for the function, specifying the class Ptype.
      \param m is a py::module object created by PYBIND11_MODULE.
      \param PT_str is a string used to create a unique Python class-name based on the Ptype.

      \return void

      \sa UserParticleBoundaryFunctions.h, UserParticleBoundaryFunctions.cpp

     */
    template <Ptype PT>
    void makeParticleMeshBoundaryConditions(py::module &m, std::string const & PT_str)
    {
      using PMBC = ParticleMeshBoundaryConditions<PT>;
      // Make a class name with the particle structure type (Ptype) string appended to the
      // string "ParticleMeshBoundaryConditions_"
      std::string pyclass_name = std::string("ParticleMeshBoundaryConditions_") + PT_str;
      
      // Create the Python binding for this class
      py::class_<PMBC>(m, pyclass_name.c_str())
        
        // The C++ ctor args types are template parameters in py::init<>()
        // Note that if a py::arg() specification is needed (e.g., to specify a
        // default value), then every argument has to have a py::arg().
        .def(py::init<std::vector<std::string>&, py::object&, UserParticleBoundaryFunctions<PT>&, bool>(), py::arg("species_names"), py::arg("pmesh_M"), py::arg("userParticleBoundaryFunctions"), py::arg("print_flag") = false);      
      // cf. ~/workspace/dolfin/python/src/geometry.cpp for __getitem__, if needed.
        
    } // void makeParticleMeshBoundaryConditions
    
  } // namespace is anonymous

}
#endif
