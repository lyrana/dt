/*! \file particle.h

  \brief This file has the source code for a C++ implementation of particle advance.

  \namespace dnt

  template <Ptype PT>
    void move_neutral_species(py::object species_name, py::object ctrl, py::object particle_P)

  template <Ptype PT>
    void move_charged_species_in_uniform_fields(py::object species_name, py::object ctrl, py::object particle_P)

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
#include "MeshEntityArrays.h"

#include "dolfin_functions.h"

//using namespace std; 

//template<typename PS>
//using SegList = std::vector<py::array_t<PS,0>>;

namespace py = pybind11;

// Static scratch space
  
double particlePosition[3];
// A scratch array that can hold: x,y,z, x0,y0,z0 (or subset)
double pCoord2[6];
// A scratch array that can hold: dx,dy,dz (or subset)
double dx[3];
// 
dolfin::Cell* pcellPtr;

namespace dnt
{

  // BEGINDEF: void move_neutral_species(py::object species_name, py::object ctrl, py::object particle_P)
  //! Advance a neutral particle species by one time-increment on a mesh.
  /*!          

    Compute change in position in time dt. Use an explicit method to calculate the
    final position.  If a particle leaves its initial cell, the cell that the
    particle moves to is calculated by finding what facet the particle crosses, and
    looking up what the neighbor cell is. This is repeated until the cell
    containing the final position is found.

    \param ctrl is a DTcontrol_C object.
    \param sap is a SegmentedArrayPair containing the particle data for one neutral species.
    \param mesh is a dolfin Mesh that the particles move through.
    \cell_vertex_dict is a dictionary that returns the vertex ids, given a cell id.

    :cvar double dtRemaining: The time a particle has left to move in the new cell
                              the particle has entered.
    :cvar int pDim: Number of spatial coordinates in the particle location.
    :cvar double tStart: The time at which a particle starts its move in the
                                 current cell.

  */
  //  template <Ptype PT>
  template <Ptype PT, size_t N_CELL_FACETS>
    void move_neutral_species(py::object particle_P, py::object species_name, py::object ctrl)
  {

    using namespace pybind11::literals;
        
    std::cout << "Hello from move_neutral_species@0" << std::endl;
    std::cout << "This is species " << std::string(py::str(species_name)) << std::endl;

    // Get attributes needed from Particle_C arg
    auto pmesh_M = particle_P.attr("pmesh_M");
    auto mesh = pmesh_M.attr("mesh").cast<dolfin::Mesh>();
    auto tDim = mesh.topology().dim();
    auto meshEntityArrays = pmesh_M.attr("mea_object").cast<MeshEntityArrays<N_CELL_FACETS> *>();
    
    auto particleBoundaryMarker = pmesh_M.attr("particle_boundary_marker").cast<dolfin::MeshFunction<size_t>>();
    
    auto NO_CELL = pmesh_M.attr("NO_CELL").cast<int>();
    auto NO_FACET = pmesh_M.attr("NO_FACET").cast<int>();
    
    // auto traj_T = particle_P.attr("traj_T").cast<py::capsule>();
    // auto traj_T = particle_P.attr("traj_T").cast<nullptr_t>();
    // There's no cast<>() here as we only need to check if traj_T is None.
    auto traj_T = particle_P.attr("traj_T");
    // Get parameters needed from ctrl arg
    auto dt = ctrl.attr("dt").cast<double>(); // This does a COPY.
    auto step = ctrl.attr("timeloop_count").cast<int>();
    auto time = ctrl.attr("time").cast<double>();
    auto MAX_FACET_CROSS_COUNT = ctrl.attr("MAX_FACET_CROSS_COUNT").cast<size_t>();
    
    //    std::cout << "NO_FACET " << NO_FACET << std::endl;
    
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

    //    std::cout << "pDim " << pDim << std::endl;

    //1 auto pseg_arr = particle_P.attr("pseg_arr").cast<py::dict>();

    //    std::cout << "Hello from move_neutral_species@1" << std::endl;
    
    // Get the data for the species to be advanced
    auto pseg_arr_map = particle_P.attr("pseg_arr").cast<std::map<std::string, SegmentedArrayPair<PT> *>>();
    auto sap = pseg_arr_map[std::string(py::str(species_name))];

    // Access the SAP without constructing a std::map. Note that sap is a pointer.
    //1 auto sap = pseg_arr[species_name].cast<SegmentedArrayPair<PT> *>();

    const py::ssize_t segmentLength = sap->get_segment_length();
    //    auto mesh_coordinates = mesh.coordinates();
    
    // Start a loop over the sap segments. segTuple contains (npSeg, psegIn, psegOut)
    py::tuple segTuple = sap->init_inout_loop(true); // 'true' returns the data pointers for
                                                     // the first segment pair, instead of the
                                                     // Numpy arrays.
    auto npSeg = segTuple[0].cast<py::ssize_t>();

    // Get the particle array pointers from the returned py::capsule's
    // psegIn:
    auto psegInCap = segTuple[1].cast<py::capsule>();
    std::cout << "The name of the capsule is " << psegInCap.name() << std::endl;
    Pstruct<PT>* psegIn = psegInCap; // Declares the type of psegIn
    std::cout << "psegIn[0].x_ = " << psegIn[0].x_ << std::endl;
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

    std::cout << "Hello from move_neutral_species@2" << std::endl;
    
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
                py::tuple segTuple = sap->get_next_out_segment(true);                
                // auto psegOutCap = segTuple[0].cast<py::capsule>();
                // Pstruct<PT>* psegOut = psegInCap;
                psegOut = segTuple[0].cast<py::capsule>(); // This compresses the above 2 lines into 1
                
                // Using Numpy structured arrays instead:
                /*
                py::tuple segTuple = sap->get_next_out_segment();
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
            //            std::cout << "ipOut " << ipOut << std::endl;
            
            psegOut[ipOut].x0_ = psegOut[ipOut].x_;
            psegOut[ipOut].y0_ = psegOut[ipOut].y_;
            psegOut[ipOut].z0_ = psegOut[ipOut].z_;
            // 2. Move the particle. NON-RELATIVISTIC: v and u are the same
            psegOut[ipOut].x_ += psegOut[ipOut].ux_*dt;
            psegOut[ipOut].y_ += psegOut[ipOut].uy_*dt;
            psegOut[ipOut].z_ += psegOut[ipOut].uz_*dt;

            auto pCellIndex = psegOut[ipOut].cell_index_;

            std::cout << "move_neutral_species. step: " << step << " ipOut: " << ipOut << " x_: " << psegOut[ipOut].x_ << std::endl;
            
            // Loop until the particle is in the current cell
            auto mLastFacet = NO_FACET;
            size_t facetCrossCount = 0;
            auto tStart = time - dt;
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
            // dolfin::Cell pcell(mesh, pCellIndex);
            // dolfin::Cell* pcellPtr = new dolfin::Cell(mesh, pCellIndex);
            pcellPtr = new dolfin::Cell(mesh, pCellIndex);
            //            auto vertices = (*pcellPtr).entities(0);
            auto vertices = pcellPtr->entities(0);

            // std::cout << "vertices " << vertices[0] << ", " << vertices[1] << std::endl;
            
            particlePosition[0] = psegOut[ipOut].x_;
            particlePosition[1] = psegOut[ipOut].y_;
            particlePosition[2] = psegOut[ipOut].z_;
            //            std::cout << "Calling is_inside_vertices" << std::endl;
            
            while (! is_inside_vertices(mesh, vertices, particlePosition))
              {
                std::cout << "Hello from call to is_inside_vertices" << std::endl;
                /*
                  The particle has left this cell.  We need to track it across each
                  facet in case there's a boundary-condition on that facet.
                */
                facetCrossCount += 1;
                // Check for an abnormal number of facet crossings:
                if (facetCrossCount > MAX_FACET_CROSS_COUNT)
                  {
                    std::cout << "DnT:particle.h:move_neutral_species: !!! MAX_FACET_CROSS_COUNT exceeded!!!" << std::endl;
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
                // auto facetNormalVector = facetTupl[2].cast<std::array<double, 3>>(); // a 3-vector
                auto facetNormalVector = facetTupl[2].cast<py::array_t<double>>();
                
                //                size_t cellFNVdim = 3; // Number of doubles per facet-normal vector (FNV).
                //                auto pyFacetNormalVector = py::array_t<double>(cellFNVdim, facetNormalVector); // Create a Numpy array and give it values.
 
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
                    //                    int mFacet = mFacetArray[cFacet]; // The mesh-level index of the facet just crossed.
                    auto mFacet = mFacetArray[cFacet]; // The mesh-level index of the facet just crossed.

                    // mFacet should never be the same as the last facet crossed: check this
                    if (mFacet == (size_t)mLastFacet) // If the particle has crossed the same facet twice in succession, there's an error:
                      {
                        std::string errorMsg = "In move_neutral_species(): The mesh index of the facet crossed is "  + std::to_string(mFacet) + ", the same as the last facet crossed. This should not happen since the particle cannot cross the same facet twice in one move!";
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
                        if ((psegOut[ipOut].bitflags_ & Pstruct<PT>::TRAJECTORY_FLAG) != 0b0) // If this is a trajectory particle.
                          {
                            // Get the storage index that currently identifies this
                            // particle in the trajectory list of particles.
                            auto fullIndex = sap->get_full_index(ipIn, "in");
                            //return py::dict("x"_a=x_, "x0"_a=x0_, "ux"_a=ux_, "weight"_a=weight_, "bitflags"_a=bitflags_, "cell_index"_a=cell_index_, "unique_ID"_a=unique_ID_, "crossings"_a=crossings_);
                            //bool facet_crossing;
                            particle_P.attr("record_trajectory_datum")(species_name, ipOut, fullIndex, step, tStart, "facet_crossing"_a=true); // Does this work?
                    /*         self.record_trajectory_datum(species_name, psegOut[ipOut], fullIndex, step, tStart, facet_crossing=True); */
                          }
                        // A reference to dx[] is available in the BC function class.
                        // double dx_fraction;
                        //particle_P.attr("pmesh_bcs").attr("bc_function_dict")[facValue][species_name](psegOut[ipOut], species_name, mFacet, dx_fraction=dxFraction, facet_normal=facetNormal);
                        // HERE2 Need to cast this first.
                        // ?? Test for is_none() first?
                        auto bcFunctionDict = particle_P.attr("pmesh_bcs").attr("bc_function_dict").cast<py::dict>();
                        //doesn't work: auto bcFunction = bcFunctionDict[std::to_string(facValue)].cast<py::dict>();
                        // A py::dict needs a py::str key value. Convert facValue to a string first.
                        auto bcFunction = bcFunctionDict[py::str(std::to_string(facValue))].cast<py::dict>();
                        // A py::dict needs a py::str key value
                        // bcFunction[py::str(species_name)](ipOut, species_name, mFacet, "dx_fraction"_a=dxFraction, "facet_normal"_a=pyFacetNormalVector);
                        bcFunction[py::str(species_name)](ipOut, species_name, mFacet, "dx_fraction"_a=dxFraction, "facet_normal"_a=facetNormalVector);                        
                    /*     self.pmesh_bcs.bc_function_dict[facValue][species_name](psegOut[ipOut], species_name, mFacet, dx_fraction=dxFraction, facet_normal=facetNormal); */
                        // particle_P.attr("pmesh_bcs").attr("bc_function_dict")[facValue][species_name](ipOut, species_name, mFacet, "dx_fraction"_a=dxFraction, "facet_normal"_a=facetNormalVector);
                      }

                    // Look up the cell index of the new cell.
                    auto pCellIndexNew = meshEntityArrays->get_cell_neighbors(pCellIndex)[cFacet];

                    /* If the particle has left the mesh, and has been deleted, end
                       the search.  If it hasn't been deleted (e.g., reflected),
                       continue tracking it.
                    */
                    if (pCellIndexNew == NO_CELL)
                      {
                        if ((psegOut[ipOut].bitflags_ & Pstruct<PT>::DELETE_FLAG) == 0b1)
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
                //                vertices = (*pcellPtr).entities(0);
                vertices = pcellPtr->entities(0);
                
                particlePosition[0] = psegOut[ipOut].x_;
                particlePosition[1] = psegOut[ipOut].y_;
                particlePosition[2] = psegOut[ipOut].z_;
              } // END: while (! is_inside_vertices(mesh, vertices, particlePosition))
            delete pcellPtr; // Not needed past this point.
//HERE2
// Add back missing code: // Record the number of facet-crossings
                          // Check for deletion
                          // Treat trajectory particles

            // Record the number of facet-crossings
            psegOut[ipOut].crossings_ = facetCrossCount;

            // For testing only:
            particle_P.attr("update_trajectory_particleId")(species_name, ipIn, ipOut);
            
            // Check that this particle has not been deleted before incrementing the
            // "out" particle counter.

            // test-only
            std::cout << "move_neutral_species2 step: " << step << " ipOut: " << ipOut << " x_: " << psegOut[ipOut].x_ << std::endl;            
            particle_P.attr("remove_trajectory_particleId")(species_name, ipIn, ipOut, step, time, dt);

            
            if ((psegOut[ipOut].bitflags_ & Pstruct<PT>::DELETE_FLAG) == 0b0) // If this particle has not been deleted...
              {
                /*
                  If particle indices in the "out" SA have changed (i.e., indexChange
                  is 'true') due to deletions, then update this particle's index where
                  needed.  E.g., if this is a trajectory particle, update its SA
                  index in the list of trajectory particles.
                */
                if (!traj_T.is_none())
                  {
                    if ((psegOut[ipOut].bitflags_ & Pstruct<PT>::TRAJECTORY_FLAG) != 0b0) // If this is a trajectory particle.
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
                    if ((psegIn[ipIn].bitflags_ & Pstruct<PT>::TRAJECTORY_FLAG) != 0b0)
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
            psegIn = segTuple[1].cast<py::capsule>(); // This is a compressed expression: there's a cast performed to make psegIn, which has it's type from above (Pstruct<PT>* psegIn)
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
    std::cout << "particleCount: " << particleCount << std::endl;
    sap->set_number_of_items("out", particleCount);
  }
  // ENDDEF: void move_neutral_species(py::object species_name, py::object ctrl, py::object particle_P)
  
  // BEGINDEF: void move_charged_species_in_uniform_fields(py::object species_name, py::object ctrl, py::object particle_P)
  //! Advance a charged-particle species one timestep in uniform fields.
  /*!          
    Apply the electric field ctrl.E0 to particles of the given species for a
    time interval ctrl.dt.  Compute the resulting change in particle
    velocities and positions.

    NB: There is no implementation of boundary-conditions in this function, so
    particles do not get deleted.

    \param species_name is the string name of the charged-particle species.
    \param ctrl is a DTcontrol_C object.
    \param sap is a SegmentedArrayPair containing the particle data for one neutral species.

  */
  template <Ptype PT>
  //  void move_charged_species_in_uniform_fields(py::str species_name, py::object ctrl, py::object particle_P)
  void move_charged_species_in_uniform_fields(py::object species_name, py::object ctrl, py::object particle_P)
  {

    //    std::cout << "Hello from move_charged_species_in_uniform_fields@0" << std::endl;

    // Get parameters needed from ctrl arg
    auto dt = ctrl.attr("dt").cast<double>(); // This does a COPY.
    auto E0x = ctrl.attr("E0").attr("x").cast<double>();
    auto E0y = ctrl.attr("E0").attr("y").cast<double>();
    auto E0z = ctrl.attr("E0").attr("z").cast<double>();

    // Get attributes needed from Particle_C arg
    auto qom = particle_P.attr("qom")[species_name].cast<double>();
    auto qmdt = qom*dt;
          
    auto pDim = particle_P.attr("particle_dimension").cast<int>();

    //    std::cout << "pDim " << pDim << std::endl;

    //    std::cout << "Hello from move_charged_species_in_uniform_fields@1" << std::endl;

    // Get the data for the species to be advanced
    //    auto pseg_arr_map = particle_P.attr("pseg_arr").cast<std::map<std::string, SegmentedArrayPair<PT> *>>();
    //    auto sap = pseg_arr_map[std::string(py::str(species_name))];
    // This is equivalent to the above two lines:
    // (pseg_arr[species_name] is a pointer to the SAP?)
    auto sap = particle_P.attr("pseg_arr")[species_name].cast<SegmentedArrayPair<PT> *>();
    
    // Access the SAP without constructing a std::map. Note that sap is a pointer.
    //1 auto sap = pseg_arr[species_name].cast<SegmentedArrayPair<PT> *>();

    const py::ssize_t segmentLength = sap->get_segment_length();
    //    auto mesh_coordinates = mesh.coordinates();
    
    // Start a loop over the sap segments. segTuple contains (npSeg, psegIn, psegOut)
    py::tuple segTuple = sap->init_inout_loop(true); // 'true' returns the data pointers for
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

    //    std::cout << "Hello from move_charged_species_in_uniform_fields@2" << std::endl;
    
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
                py::tuple segTuple = sap->get_next_out_segment(true);                
                // auto psegOutCap = segTuple[0].cast<py::capsule>();
                // Pstruct<PT>* psegOut = psegInCap;
                psegOut = segTuple[0].cast<py::capsule>(); // This compresses the above 2 lines into 1
                
                // Using Numpy structured arrays instead:
                /*
                py::tuple segTuple = sap->get_next_out_segment();
                auto psegOutNumpyArray = segTuple[0].cast<py::array_t<Pstruct<PT>,0>>();
                const auto psegOutNumpyArray_info = psegOutNumpyArray.request();
                const auto psegOut = static_cast<Pstruct<PT>*>(psegOutNumpyArray_info.ptr);
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
    sap->set_number_of_items("out", particleCount);
  }
  // ENDDEF: void move_charged_species_in_uniform_fields(py::object species_name, py::object ctrl, py::object particle_P)
  
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

  void move_neutral_species(SegmentedArrayPair<Ptype::cartesian_xyz>& sap)
  {
    
  }
  */
  
} // namespace dnt

#endif
