# UserParticles Module

"""
These are static (like global) functions. No SELF variable!
"""
__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
__all__ = ['UserParticleDistributions_C',
           'UserParticleBoundaryFunctions_C',
           'UserParticleSourceFunctions_C',
           ]

import sys
import math
import numpy as np_M

from Dolfin_Module import Mesh_C
from Particle_Module import Particle_C
from Particle_Module import ParticleSpecies_C
#from Particle_Module import ParticleMeshBoundaryConditions_C
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

### The functions below are needed for particles that exist at the start of the
### calculation.

class UserParticleDistributions_C(object):
    """UserParticleDistributions_C is to be edited by the user to specify initial particle
       distributions and other attributes.  The units are MKS by default (i.e., if no
       conversion factor is applied), but any units can be used provided the conversion to
       MKS is available in the UserUnits_M.py module.

       .. note:: Position-coordinate components and force components are set using the
                 *position_coordinates* and *force_components* attributes, respectively,
                 of the ParticleInput_C object.  These are usually specified by the user
                 in the main source file.
    """

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

        ## First particle

        # Particle moves in -x direction:
        (x0,)  = (9.5,)
        (ux0,) = (-2.0,)

        weight0 = 2.0 # number of ions per macroparticle
        bitflags0 = 0b0 # bit flags variable
        # Turn ON trajectory flag.
        bitflags0 = bitflags0 | Particle_C.TRAJECTORY_FLAG
        cell_index = Mesh_C.NO_CELL
        unique_ID = Particle_C.UNIQUE_ID_COUNTER; Particle_C.UNIQUE_ID_COUNTER += 1
        
        # Trim the number of coordinates here if needed to match "position_coordinates" variable in ParticleInput_C
        p0 = (x0, x0, ux0, weight0, bitflags0, cell_index, unique_ID)

        ## Second particle

        # Particle moves in +x direction:
        (x1,) = (-9.5,)
        (ux1,) = (2.0,)
        weight1 = 3.0
        bitflags1 = 0b0 # bit flags variable
        # Turn ON trajectory flag.
        bitflags1 = bitflags1 | Particle_C.TRAJECTORY_FLAG 
        unique_ID = Particle_C.UNIQUE_ID_COUNTER; Particle_C.UNIQUE_ID_COUNTER += 1
        
        # Trim the number of coordinates here to match "position_coordinates" variable in ParticleInput_C
        p1 = (x1, x1, ux1, weight1, bitflags1, cell_index, unique_ID)

#        particle_list = (p0,) # Need the final comma , for one particle.
#        particle_list = (p1,) # Need the final comma , for one particle.
        particle_list = (p0, p1)

        np = len(particle_list)

        return np, particle_list
#    def neutral_H(type):ENDDEF

#class UserParticleDistributions_C(object):ENDCLASS

class UserParticleBoundaryFunctions_C(object):
    """UserParticleBoundaryFunctions_C implements callback functions
       (boundary conditions) for kinetic particles crossing marked
       mesh facets.

       See Particle_Module::ParticleMeshBoundaryConditions_C for naming
       scheme.
    """

    @staticmethod
    def default_bc(p, speciesName, facetIndex):
        """Global default boundary condition for all species.

           :param p: the record of the particle that crossed xmin.
           :param speciesName: the species that particle p belongs to.
           :param facetIndex: the facet crossed by particle p.
        """
        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print "Called", fncName
        print "Particle is", p, "at facet", facetIndex

        # Set the delete flag
        p['bitflags'] = p['bitflags'] | Particle_C.DELETE_FLAG

        # Count the number/charge/energy of deleted particles

        return
#    def default_bc(p, speciesName, facetIndex):ENDDEF

    # For Cartesian coordinates

#class UserParticleBoundaryFunctions_C(object):
    @staticmethod
    def default_bc_at_xmin(p, speciesName, facetIndex):
        """Default boundary condition for particles incident on xmin.

           :param p: the record of the particle that crossed xmin.
           :param speciesName: the species that particle p belongs to.
           :param facetIndex: the facet crossed by particle p.
        """
        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():'
        print fncName, "invoked by particle", p, "of species", speciesName

        # Set the delete flag
        p['bitflags'] = p['bitflags'] | Particle_C.DELETE_FLAG

        # Count the number/charge/energy of deleted particles

        return
#    def default_bc_at_xmin(p, speciesName, facetIndex):ENDDEF
    
#class UserParticleBoundaryFunctions_C(object):
    @staticmethod
    def default_bc_at_xmax(p, speciesName, facetIndex):
        """Default boundary condition for particles incident on xmax.

           :param p: the record of the particle that crossed xmin.
           :param speciesName: the species that particle p belongs to.
           :param facetIndex: the facet crossed by particle p.
        """
        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print "Called", fncName

        return
#    def default_bc_at_xmax(p, speciesName, facetIndex):ENDDEF
    
#class UserParticleBoundaryFunctions_C(object):
    @staticmethod
    def bc_at_xmin_for_neutral_H(p, speciesName, facetIndex):
        """Boundary condition for neutral_H incident on xmin.

           :param p: the record of the particle that crossed xmin.
           :param speciesName: the species that particle p belongs to.
           :param facetIndex: the facet crossed by particle p.
        """
        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():'
        print fncName, "invoked by particle", p, "of species", speciesName

        # Set the delete flag
        p['bitflags'] = p['bitflags'] | Particle_C.DELETE_FLAG


    # For spherical coordinates

#class UserParticleBoundaryFunctions_C(object):
    @staticmethod
    def default_bc_at_rmin(p, speciesName, facetIndex):
        """Default boundary condition for particles incident on rmin.

           :param p: the record of the particle that crossed rmin.
           :param speciesName: the species that particle p belongs to.
           :param facetIndex: the facet crossed by particle p.
        """
        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():'
#        print fncName, "invoked by particle", p, "of species", speciesName

        # Set the delete bit
        p['bitflags'] = p['bitflags'] | Particle_C.DELETE_FLAG

        # Count the number/charge/energy of deleted particles

        return
#    def default_bc_at_rmin(p, speciesName, facetIndex):ENDDEF
    
#class UserParticleBoundaryFunctions_C(object):
    @staticmethod
    def default_bc_at_rmax(p, speciesName, facetIndex):
        """Default boundary condition for particles incident on rmax.

           :param p: the record of the particle that crossed rmin.
           :param speciesName: the species that particle p belongs to.
           :param facetIndex: the facet crossed by particle p.
        """
        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
#        print "Called", fncName

        # Set the delete flag
        p['bitflags'] = p['bitflags'] | Particle_C.DELETE_FLAG

        return
#    def default_bc_at_rmax(p, speciesName, facetIndex):ENDDEF
    
#class UserParticleBoundaryFunctions_C(object):
    @staticmethod
    def bc_at_rmin_for_neutral_H(p, speciesName, facetIndex):
        """Boundary condition for neutral_H incident on rmin.

           :param p: the record of the particle that crossed rmin.
           :param speciesName: the species that particle p belongs to.
           :param facetIndex: the facet crossed by particle p.
        """
        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():'
        print fncName, "invoked by particle", p, "of species", speciesName

        # Set the delete flag
        p['bitflags'] = p['bitflags'] | Particle_C.DELETE_FLAG

        # Count the number/charge/energy of deleted particles

        return
#    def bc_at_rmin_for_neutral_H(p, speciesName, facetIndex):ENDDEF

#class UserParticleBoundaryFunctions_C(object): ENDCLASS

class UserParticleSourceFunctions_C(object):
    """UserParticleSourceFunctions_C implements particle sources

       See Particle_Module::ParticleMeshBoundaryConditions_C for naming
       scheme.
    """

## XX these are supposed to be SOURCE functions.
    
#class UserParticleSourceFunctions_C(object):
    @staticmethod
    def default_bc(p, speciesName, facetIndex):
        """Global default boundary condition for all species.

           :param p: the record of the particle that crossed xmin.
           :param speciesName: the species that particle p belongs to.
           :param facetIndex: the facet crossed by particle p.
        """
        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print "Called", fncName

        return

#class UserParticleSourceFunctions_C(object):
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

class SpecialElectrons_C(ParticleSpecies_C):
    """SpecialElectrons_C implements a source of electrons.

       See Particle_Module::ParticleMeshBoundaryConditions_C for naming
       scheme.
    """

    def __init__(self, charge, mass, dynamics):
        """Initialize an SpecialElectrons_C instance.

           :param float charge: Physical charge of one particle of the species
           :param float mass: Physical mass of one particle of the species
           :param str dynamics: Determines the type of algorithm used to push
                                the particles.  Can be either 'explicit' or 'implicit'.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

        # Pass the map_tol argument on to the base class
        super(self.__class__, self).__init__(charge, mass, dynamics)

        return
#    def __init__(self, charge, mass, dynamics):ENDDEF

#class SpecialElectrons_C(ParticleSpecies_C):
    def add_particles(self, cell_list):
        """Adds electrons to a subdomain defined by a cell list.

           :param cell_list: The cells belonging to the mesh subdomain.

           :returns: None

       """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

        return
#    def add_particles(self, cell_list):ENDDEF

#class SpecialElectrons_C(object):ENDCLASS
