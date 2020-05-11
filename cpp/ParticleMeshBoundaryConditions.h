/*! \file ParticleMeshBoundaryConditions.h

  \brief This file defines the C++ class ParticleMeshBoundaryConditions, which contains
  the map from mesh-facet tags to the call-back functions treating boundary-crossing
  particles.

  This is the source code for a C++ implementation of the Python class
  ParticleMeshBoundaryConditions_C.

  \namespace dnt
  \sa user_particle_boundary_functions_solib.cpp UserParticleBoundaryFunctions.h, particle.h

*/

#ifndef PARTICLEMESHBOUNDARYCONDITIONS_H
#define PARTICLEMESHBOUNDARYCONDITIONS_H

#include <iostream>
// uncomment to disable assert()
// #define NDEBUG
#include <cassert>
#include "Pstruct.h"

namespace dnt
{
  /*! \class ParticleMeshBoundaryConditions
    \brief This class contains a map from mesh-facet tags to callback functions.

    This is a C++ version of the Python class ParticleMeshBoundaryConditions_C.
  */
  template <Ptype PT>    
    class ParticleMeshBoundaryConditions
    {
    private:
      // This map (dictionary) contains a pointer to boundary callback
      // function for each [boundary tag][species] pair (i.e., the key is a 2D object).
      std::map<std::pair<int, std::string>, void (*)()> bc_function_dict;
      
    public:
      //! The one and only user-defined ParticleMeshBoundaryConditions ctor.
      /*!

        Create a map from facet tags to callback functions. This is how particle boundary
        conditions are implemented. The facet tags are the values of a dolfin FacetFunction
        defined on mesh facets.
      
        The following function naming-scheme is used for the callback functions:
          default_bc(): The default function for all species.
          default_bc_at_name(): The default called for all species at
          the boundary 'name'.
          bc_at_name_for_species: The function called for 'species'
          crossing 'name'.

        The most specific function found for a given boundary and
        species is used.

        \param species_names: Species names in an std::vector<std::string>

        \param pmesh_M: The Mesh_C object use by the particle advance.

        \param userPBFs: A UserParticlesBoundaryFunctions object that contains the callback
                         functions for particle boundary-conditions.


        \sa Ptype

      */
      ParticleMeshBoundaryConditions(std::vector<std::string>& species_names, py::object pmesh_M, UserParticleBoundaryFunctions<PT> userPBFs, bool print_flag = False)
        {

          // Set local names from passed parameters
          // pBDict = pmesh_M.particle_boundary_dict
          auto pBDict = particle_P.attr("particle_boundary_dict").cast<std::map<std::string, int>>();
        /*
          We want to associate the names of the user-supplied boundary functions with the
          int values of the facet tags.  First, we swap the BC dictionary keys and values:
          e.g., 'xmin': 1 to 1: 'xmin'. This lets us use the value of the facet function,
          a uint (size_t), as an index to get the string name of the boundary.
        */
        // pBDictInv = {v: k for k, v in pBDict.items()}
          std::map<int, std::string> pBDictInv; // The reversed dictionary
          for (std::map<int, std::string>::iterator it = pBDict.begin(); it != pBDict.end(); ++it)
            pBDictInv[it->second] = it->first;
 
          // Get a list of the key values in pBDictInv:
          // particleBoundaryTags = list(pBDictInv.keys())
          std::vector<int> particleBoundaryTags;
          for(std::map<int, std::string>::iterator it = pBDictInv.begin(); it != pBDictInv.end(); ++it)
            {
              particleBoundaryTags.push_back(it->first);
              std::cout << it->first << std::endl;
            }        

          // Initialize the dictionary that will contain the name of a boundary callback
          // function for each [boundary tag][species] pair (i.e., the key is a 2D object).
          // self.bc_function_dict = dict((intTag, dict((sp, None) for sp in speciesNames)) for intTag in particleBoundaryTags)
          // Use a std::pair as the key [boundary-tag][species-name] index:

          //std::vector<void (*)()> vectoroffunctions;
 
          //std::map<std::pair<int,int>, int> myMap;
          //myMap[std::make_pair(10,20)] = 25;
          //std::cout << myMap[std::make_pair(10,20)] << std::endl;
 
          // Find the default boundary function for all boundaries and all species (if there
          // is one), in the UserParticleBoundaryFunctions object passed in.
          // These local variables are used in the search:
          std::string bcFunctionName;
          void (*)() bcFunctionPtr;
          void (*)() bcGlobalDefaultFunctionPtr;
          void (*)() bcBoundaryDefaultFunctionPtr;
          // Create a short alias for the BC function map passed in:
          std::map<std::string, void (*)()>& bcFunctionMap = userPBFs<PT>.bc_function_map;
 
          // Check if a global default BC function has been defined:
          bcFunctionName = "default_bc";
          std::map<std::string, void (*)()>::iterator it = bcFunctionMap.find(bcFunctionName);
          if(it != bcFunctionMap.end())
            {
              bcGlobalDefaultFunctionPtr = it->second;
            }
          else
            {
              bcGlobalDefaultFunctionPtr = nullptr;
            }
 
          // Now loop on particle boundary tags and on particle species names to find the
          // most specific BC. Overwrite the default function with a function specific
          // to each boundary, if there is one.
          for (std::vector<int>::iterator intTag = particleBoundaryTags.begin(); intTag != particleBoundaryTags.end(); intTag++)
          {
            bcFunctionName = std::string("default_bc_at_") + pBDictInv[intTag];
            // Look for a callback function specific to this boundary. If one is found,
            // overwrite the global default callback.
            std::map<std::string, void (*)()>::iterator it = bcFunctionMap.find(bcFunctionName);
            if(it != bcFunctionMap.end())
              {
                bcBoundaryDefaultFunctionPtr = it->second;
              }
            else
              {
                bcBoundaryDefaultFunctionPtr = bcGlobalDefaultFunctionPtr;
              }

            for (std::vector<std::string>::iterator sp = species_names.begin(); sp != species_names.end(); sp++)
              {
                bcFunctionName = std:string("bc_at_") + pBDictInv[intTag] + "_for_" + sp;
                // Look for a function specific to this boundary and species. If there is
                // one, then overwrite the default for this boundary.
                std::map<std::string, void (*)()>::iterator it = bcFunctionMap.find(bcFunctionName);
                if(it != bcFunctionMap.end())
                  {
                    bcFunctionPtr = it->second;
                  }
                else
                  {
                    bcFunctionPtr = bcBoundaryDefaultFunctionPtr;
                  }
                if (bcFunctionPtr == nullptr)
                  {
                    // There's no boundary-condition for this species at this boundary.
                    std::cout << "{ParticleMeshBoundaryConditions.h}ParticleMeshBoundaryConditions(): No callback function was specified for " << pBDictInv[intTag] << "/" << sp << std::endl;
                    std::exit(EXIT_FAILURE);
                  }
                else if (print_flag == true)
                  {
                    std::cout << "{ParticleMeshBoundaryConditions.h}ParticleMeshBoundaryConditions(): Callback function for " << pBDictInv[intTag] << "/" << sp << " is " << bcFunctionName << std::endl;
                  }
                // Finally, store the most-specific callback in the dictionary:
                bc_function_dict[std::make_pair(intTag,sp)] = bcFunctionPtr;
              } // Loop on species
          } // Loop on boundary tags
          
        } // ParticleMeshBoundaryConditions(std::vector<std::string>& species_names, py::object pmesh_M, UserParticleBoundaryFunctions<PT> userPBFs, bool print_flag = False): ENDDEF

    } // class ParticleMeshBoundaryConditions: ENDCLASS
  
} // namespace dnt
        
#endif
