# UserParticles Module

"""
These are static (like global) functions. No SELF variable!
"""
__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
__all__ = ['UserParticleDistributions_C.drifting_Hplus',
           ]

import sys
import math
import numpy as np_m

from Dolfin_Module import Mesh_C
from Particle_Module import Particle_C
import UserUnits_Module as U_M

#
# these may need info like the mesh elements, etc.  Could get this by
# importing the mesh?  If the distribution is single particles, just
# need to add particles to their elements and to the storage array.
#
# Maybe they just have the args they need
#

# Select the unit system to be used for input parameters.
Convert = U_M.MyPlasmaUnits_C

#STARTCLASS
class UserParticleDistributions_C(object):
    """UserParticleDistributions_C is to be edited by the user to
       specify initial particle distributions.  The units are MKS by
       default (i.e., if no conversion factor is applied), but any
       units can be used provided the conversion to MKS is available
       in the UserUnits_M.py module.
    """

    # The spatial particle coordinates to be saved in particle storage.
#    position_coordinates = ['x', 'y', 'z']

# Initialize a few particles from each species
# 
    @staticmethod
    def neutral_H(type):
        """Neutral hydrogen drifting at constant speed.
           
           See ParticleNoFields.ods for test results.
        """
        # Check that the caller knows this function returns 'listed' particles:
        if type != 'listed':
            error_msg = "neutral_H() called with type '" + str(type) + "' but it is of type 'listed'"
            sys.exit(error_msg)

        # Set initial phase-space coordinates and macroparticle weight for each particle

        # First particle
        # Particle moves in -x direction at constant y, z:
        (x0, y0, z0) = (9.5, -9.5, 0.0)
        (ux0, uy0, uz0) = (-2.0, 0.0, 0.0)

        weight0 = 2.0 # number of ions per macroparticle
        bitflags0 = 0b0 # bit flags variable
        # Turn ON trajectory flag.
        bitflags0 = bitflags0 | Particle_C.TRAJECTORY_FLAG
        cell_index = Mesh_C.NO_CELL
#        cell_index = -1
        unique_ID = Particle_C.UNIQUE_ID_COUNTER; Particle_C.UNIQUE_ID_COUNTER += 1
        crossings = 0

        # Trim the number of coordinates here if needed to match "position_coordinates" variable in ParticleInput_C
        p0 = (x0,y0,z0, x0,y0,z0, ux0,uy0,uz0, weight0, bitflags0, cell_index, unique_ID, crossings)

        # Second particle
        # Particle moves diagonally across the mesh:
        (x1, y1, z1) = (9.5, 9.5, 9.5)
        (ux1, uy1, uz1) = (-2.0, -2.0, -2.0)
        weight1 = 3.0
        bitflags1 = 0b0 # bit flags variable
        # Turn ON trajectory flag.
        bitflags1 = bitflags1 | Particle_C.TRAJECTORY_FLAG 
        unique_ID = Particle_C.UNIQUE_ID_COUNTER; Particle_C.UNIQUE_ID_COUNTER += 1
        crossings = 0
        
        # Trim the number of coordinates here to match "position_coordinates" variable in ParticleInput_C
        p1 = (x1,y1,z1, x1,y1,z1, ux1,uy1,uz1, weight1, bitflags1, cell_index, unique_ID, crossings)

#        particle_list = (p0,) # Need the final comma , for one particle.
#        particle_list = (p1,) # Need the final comma , for one particle.
        particle_list = (p0, p1)

        np = len(particle_list)

        return np, particle_list
#    def neutral_H(type):ENDDEF

#class UserParticleDistributions_C(object):ENDCLASS

# Does this need to descend from ParticleBoundaryConditions_C?  It just has static functions
#class UserParticleMeshBoundaryConditions_C(ParticleBoundaryConditions_C):
#class UserParticleMeshBoundaryConditions_C(object):
#STARTCLASS
class UserParticleBoundaryFunctions_C(object):
    """UserParticleBoundaryFunctions_C implements callback functions
       (boundary conditions) for kinetic particles crossing marked
       mesh facets.

       See Particle_Module::ParticleMeshBoundaryConditions_C for naming
       scheme.
    """

    @staticmethod
    def default_bc(p, speciesName, facetIndex, dx_fraction=None, facet_normal=None):
        """Global default boundary condition for all species.
        """
        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print("Called", fncName)

        return

    @staticmethod
    def default_bc_at_xmin(p, speciesName, facetIndex, dx_fraction=None, facet_normal=None):
        """Default boundary condition for particles incident on xmin.

           :param p: the record of the particle that crossed xmin.
           :param speciesName: the species that particle p belongs to.
           :param facetIndex: the facet crossed by particle p.
        """

        printInfoInvoked = False

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():'

        if printInfoInvoked is True:
            print(fncName, "invoked by particle", p, "of species", speciesName)

        # Set the delete flag
        p['bitflags'] = p['bitflags'] | Particle_C.DELETE_FLAG

        # Count the number/charge/energy of deleted particles

        return
#    def default_bc_at_xmin(p, facetIndex):ENDDEF
    
    @staticmethod
    def default_bc_at_xmax(p, speciesName, facetIndex, dx_fraction=None, facet_normal=None):
        """Default boundary condition for particles incident on xmax.
        """
        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print("Called", fncName)

        return
    
    @staticmethod
    def default_bc_at_ymin(p, speciesName, facetIndex, dx_fraction=None, facet_normal=None):
        """Default boundary condition for particles incident on ymin.
        """

        printInfoInvoked = False
        
        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():'

        if printInfoInvoked is True:        
            print(fncName, "invoked by particle", p, "of species", speciesName)

        # Set the delete flag
        p['bitflags'] = p['bitflags'] | Particle_C.DELETE_FLAG

        # Count the number/charge/energy of deleted particles

        return
    
    @staticmethod
    def default_bc_at_ymax(p, speciesName, facetIndex, dx_fraction=None, facet_normal=None):
        """Default boundary condition for particles incident on ymax.
        """
        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print("Called", fncName)

        return
    
    @staticmethod
    def bc_at_xmin_for_neutral_H(p, speciesName, facetIndex, dx_fraction=None, facet_normal=None):
        """Boundary condition for neutral_H incident on xmin.

        """

        printInfoInvoked = False
        
        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():'

        if printInfoInvoked is True:
            print(fncName, "invoked by particle", p, "of species", speciesName)

        # Set the delete flag
        p['bitflags'] = p['bitflags'] | Particle_C.DELETE_FLAG

        # Count the number/charge/energy of deleted particles

        return

#class UserParticleBoundaryFunctions_C(object): ENDCLASS

#STARTCLASS
class UserParticleSourceFunctions_C(object):
    """UserParticleSourceFunctions_C implements particle sources

       See Particle_Module::ParticleMeshBoundaryConditions_C for naming
       scheme.
    """

    @staticmethod
    def default_bc(p, speciesName, facetIndex):
        """Global default boundary condition for all species.
        """
        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print("Called", fncName)

        return

    @staticmethod
    def default_bc_at_xmin(p, speciesName, facetIndex):
#    def default_bc_at_xmin(self, p, facetIndex):
        """Default boundary condition for particles incident on xmin.

           :param p: the record of the particle that crossed xmin.
           :param speciesName: the species that particle p belongs to.
           :param facetIndex: the facet crossed by particle p.
        """
        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():'

#class UserParticleSourceFunctions_C(object):ENDCLASS
