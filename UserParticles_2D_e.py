# UserParticles Module

"""
These are static (like global) functions. No SELF variable!
"""
__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
__all__ = ['UserParticleDistributions_C.initial_electrons',
           'UserParticleDistributions_C.initial_ions',
           'UserParticleDistributions_C.few_initial_electrons',
           'UserParticleDistributions_C.few_initial_ions',
           ]

import sys
import math
import numpy as np_M

from Dolfin_Module import Mesh_C
from Particle_Module import Particle_C
from Particle_Module import ParticleMeshBoundaryConditions_C
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

class UserParticleDistributions_C(object):
    """UserParticleDistributions_C is to be edited by the user to
       specify initial particle distributions.  The units are MKS by
       default (i.e., if no conversion factor is applied), but any
       units can be used provided the conversion to MKS is available
       in the UserUnits_M.py module.
    """

    # The spatial coordinates of all particles, a class variable
#    position_coordinates = ['x', 'y',]


# Initialize a few particles from each species
# 
    @staticmethod
    def trajelectrons(type):
        """These electrons are specified as lists of coordinates and weights.
        """
        # Check that the caller knows this function returns 'listed' particles:
        if type != 'listed':
            error_msg = "trajelectrons() called with type '" + str(type) + "' but it is of type 'listed'"
            sys.exit(error_msg)

        # Set initial phase-space coordinates and macroparticle weight for each particle

        # An electron along the 45-degree line, near the outer boundary
        (x0, y0, z0) = (4.75, 0.0, 0.0)
        # Rotate through theta = pi/4 in x-y plane
# test
#        theta0 = math.pi/8.0
        theta0 = math.pi/4.0
        xnew = x0*np_M.cos(theta0)-y0*np_M.sin(theta0)
        ynew = x0*np_M.sin(theta0)+y0*np_M.cos(theta0)
        (x0, y0, z0) = (xnew, ynew, z0)
        (ux0, uy0, uz0) = (0.0, 0.0, 1000.0)
        weight0 = 1.0 # number of electrons per macroparticle
        bitflags0 = 0b00 # bit flags variable is all zeroes
        # Turn trajectory flag ON.
        bitflags0 = bitflags0 | Particle_C.TRAJECTORY_FLAG
        cell_index = Mesh_C.NO_CELL

        # Trim the number of coordinates here to match
        # "position_coordinates" variable in DTparticleInput_C
        p0 = (x0,y0, x0,y0, ux0,uy0, weight0, bitflags0, cell_index)

        # Second particle
        # An electron along the 22.5-degree line
        (x1, y1, z1) = (4.75, 0.0, 0.0)
        # Rotate through theta = pi/8 in x-y plane
# test
#        theta1 = math.pi/4.0
        theta1 = math.pi/8.0
        xnew = x1*np_M.cos(theta1)-y1*np_M.sin(theta1)
        ynew = x1*np_M.sin(theta1)+y1*np_M.cos(theta1)
        (x1, y1, z1) = (xnew, ynew, z1)
        (ux1, uy1, uz1) = (0.0, 0.0, 1000.0)
        weight1 = 2.0
        bitflags1 = 0b00 # bit flags variable is all zeroes
        # Turn trajectory flag ON.
        bitflags1 = bitflags1 | Particle_C.TRAJECTORY_FLAG
#        cell_index = -1
        cell_index = Mesh_C.NO_CELL

        # Trim the number of coordinates here to match
        # "position_coordinates" variable in DTparticleInput_C
        p1 = (x1,y1, x1,y1, ux1,uy1, weight1, bitflags1, cell_index)

        particle_list = (p0, p1)
        np = len(particle_list)

# Unit conversion                       
        return np, particle_list

    @staticmethod
    def testelectrons(type):
        """These electrons are specified as lists of coordinates and weights.
        """
        # Check that the caller knows this function returns 'listed' particles:
        if type != 'listed':
            error_msg = "testelectrons() called with type '" + str(type) + "' but it is of type 'listed'"
            sys.exit(error_msg)

        # Set initial phase-space coordinates and macroparticle weight for each particle

        # An electron at the y = 0 boundary

#        (x0, y0, z0) = (1.1, 0.1, 0.0)
#        (x0, y0, z0) = (1.0, 0.0, 0.0)
        (x0, y0, z0) = (4.75, 0.0, 0.0)
# Rotate through theta in x-y plane
        theta0 = math.pi/4.0
        xnew = x0*np_M.cos(theta0)-y0*np_M.sin(theta0)
        ynew = x0*np_M.sin(theta0)+y0*np_M.cos(theta0)
        (x0, y0, z0) = (xnew, ynew, z0)

#        print 'x0, y0, z0 =', x0, y0, z0

#        (ux0, uy0, uz0) = (3000.0, 2000.0, 1000.0)
        (ux0, uy0, uz0) = (0.0, 0.0, 1000.0)
        weight0 = 1.0 # number of electrons per macroparticle
        bitflags0 = 0b00 # bit flags variable is all zeroes
        # Turn on trajectory flag.
        bitflags0 = bitflags0 | Particle_C.TRAJECTORY_FLAG
        cell_index = Mesh_C.NO_CELL

        # Trim the number of coordinates here to match
        # "position_coordinates" variable in DTparticleInput_C
        p0 = (x0,y0, x0,y0, ux0,uy0, weight0, bitflags0, cell_index)

        # Second particle
        (x1, y1, z1) = (1.1, 0.0, 10.0)
# Rotate through theta in x-y plane
        theta1 = math.pi/4.0
        xnew = x1*np_M.cos(theta1)-y1*np_M.sin(theta1)
        ynew = x1*np_M.sin(theta1)+y1*np_M.cos(theta1)
        (x1, y1, z1) = (xnew, ynew, z1)

#        print 'x1, y1, z1 =', x1, y1, z1

#        (ux1, uy1, uz1) = (3000.0, 2000.0, 1000.0)
        (ux1, uy1, uz1) = (0.0, 0.0, 1000.0)
        weight1 = 2.0
        bitflags1 = 0b00 # bit flags variable
#        TRAJECTORY_FLAG = 0b1
        # Turn on trajectory flag.
#        bitflags1 = bitflags1 | Particle_C.TRAJECTORY_FLAG 

        # Trim the number of coordinates here to to match
        # "position_coordinates" variable in DTparticleInput_C
        p1 = (x1,y1, x1,y1, ux1,uy1, weight1, bitflags1, cell_index)

        particle_list = (p0, p1)
        np = len(particle_list)

# Unit conversion                       
#        velocity = [v*Convert.m_per_s for v in velocity]
        return np, particle_list

#class UserParticleDistributions_C(object):ENDCLASS

#class UserParticleMeshBoundaryConditions_C(ParticleBoundaryConditions_C):
class UserParticleMeshFunctions_C(object):
    """UserParticleMeshFunctions_C implements callback functions
       (boundary conditions) for kinetic particles crossing marked
       mesh facets.

       See Particle_Module::ParticleMeshBoundaryConditions_C for naming
       scheme.
    """

    @staticmethod
    def default_bc(self, p, facet_index):
        """Global default boundary condition for all species.
        """

        fncname = sys._getframe().f_code.co_name + '():'
        print "Called", fncname

        return

    @staticmethod
    def default_bc_at_rmin(self, p, facet_index):
        """Default boundary condition particles incident on rmin.
        """

        fncname = sys._getframe().f_code.co_name + '():'
        print "Called", fncname

        return
    
    @staticmethod
    def bc_at_rmin_for_testelectrons(self, p, facet_index):
        """Boundary condition for testelectrons incident on rmin.
        """

        fncname = sys._getframe().f_code.co_name + '():'
        print "Called", fncname

        return

#class UserParticleMeshFunctions_C(object): ENDCLASS
