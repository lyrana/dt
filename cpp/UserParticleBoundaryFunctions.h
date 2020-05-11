/*! \file UserParticleBoundaryFunctions.h

  \brief This file has C++ call-back functions to treat boundary-crossing particles.

  This is the source code for a C++ implementation of the Python class
  UserParticleBoundaryFunctions_C.

  \namespace dnt
  \sa user_particles_2D_e.solib.cpp

*/

#ifndef USERPARTICLEBOUNDARYFUNCTIONS_H
#define USERPARTICLEBOUNDARYFUNCTIONS_H

#include <iostream>
// uncomment to disable assert()
// #define NDEBUG
#include <cassert>
#include "Pstruct.h"

namespace dnt
{
  /*! \class UserParticleBoundaryFunctions
    \brief This class implements call-back functions that treat boundary-crossing particles.

    This is a C++ version of the Python class UserParticleBoundaryFunctions_C.

    See Particle_Module::ParticleMeshBoundaryConditions_C for the naming
    scheme for particle call-back functions.

    The class is templated on the type of the particle struct (PT) so that we can
    pass a (reference to a) particle struct as a argument to the member functions.

    \param PT is the type of particle struct that this class can operate on.
    \sa Ptype, Pstruct

  */
  template<Ptype PT>
    class UserParticleBoundaryFunctions
    {

    public:
      //! The one and only user-defined UserParticleBoundaryFunctions ctor.
      /*!

        \param position_coordinates: Example: ['x', 'y',]
        
        \sa Ptype

      */
      // Call example (sphere1D.py):
      //     userPBndFns = UserParticleBoundaryFunctions_C(particle_P.position_coordinates, particle_P.dx)
      UserParticleBoundaryFunctions(py::list& position_coordinates_arg)
        {
          for (auto item : position_coordinates_arg)
            {
              position_coordinates.push_back(item.cast<std::string>());
              // or: position_coordinates.push_back(std::string(py::str(item)));
            }
          // Create an std::map of function-names to functions
          // This has to be done manually for each callback function defined below.
          bc_function_map.insert(std::pair("default_bc", &default_bc);
 
        };
      // UserParticleBoundaryFunctions(position_coordinates, dx):ENDDEF

    private:
      std::vector<std::string> position_coordinates;
      // This is a dictionary of functions, indexed by the functions names.
      std::map<std::string, void (*)()> bc_function_map;
      
      // Scratch for manipulating the particle coordinates and velocities:
      double pcoord[3], pvel[3];
    
      //! Define the default boundary condition for all particles on all boundaries.
      /*!

        p, species_name, facet_index, dx, dx_fraction, facet_normal

        \param string in_out: Either "in" or "out" depending on whether
        we're dealing with the in or out SA.

        \return: Nothing is returned.
      */
      // Call in Particle_Module.py:
      //   self.pmesh_bcs.bc_function_dict[facValue][sn](psegOut[ipOut], sn, mFacet, dx_fraction=dxFraction, facet_normal=facetNormal)

      // Call in particle.h:
      // old: bcFunction[py::str(species_name)](ipOut, species_name, mFacet, "dx_fraction"_a=dxFraction, "facet_normal"_a=facetNormalVector);        
      //   bcFunction[species_name](ipOut, species_name, mFacet, "dx_fraction"_a=dxFraction, "facet_normal"_a=facetNormalVector);

    public:
      void default_bc(Pstruct<PT>& p, py::str& species_name, const int facet_index, const double dx[], const double dx_fraction, py::array_t<double>& facet_normal)
      {

        std::cout << "Hello from the UserParticleBoundaryFunctions::default_bc" << std::endl;

        // Set the delete flag on the particle
        p.bitflags_ = p.bitflags_ | Pstruct<PT>::DELETE_FLAG;

        // We might want to count the number/charge/energy of deleted particles here
        // before returning.
      
        return;
      }
      // void default_bc(p, species_name, facet_index, dx_fraction=None, facet_normal=None): ENDDEF

    };
  // class UserParticleBoundaryFunctions: ENDCLASS

} // namespace dnt


#endif


