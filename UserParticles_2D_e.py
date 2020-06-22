# UserParticles Module

"""
These are static (like global) functions. No SELF variable!
"""
__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
__all__ = ['UserParticleDistributions_C.trajelectrons',
           'UserParticleDistributions_C.test_electrons',
           'UserParticleDistributions_C.two_electrons',
           'UserParticleDistributions_C.few_initial_electrons',
           'UserParticleDistributions_C.few_initial_ions',

           'UserParticleBoundaryFunctions_C.default_bc',
           'UserParticleBoundaryFunctions_C.default_bc_at_rmin',
           'UserParticleBoundaryFunctions_C.bc_at_rmin_for_test_electrons',
           ]

import sys
import math
import numpy as np_m

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

#STARTCLASS
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
        xnew = x0*np_m.cos(theta0)-y0*np_m.sin(theta0)
        ynew = x0*np_m.sin(theta0)+y0*np_m.cos(theta0)
        (x0, y0, z0) = (xnew, ynew, z0)
        (ux0, uy0, uz0) = (0.0, 0.0, 1000.0)
        weight0 = 1.0 # number of electrons per macroparticle
        bitflags0 = 0b00 # bit flags variable is all zeroes
        # Turn trajectory flag ON.
        bitflags0 = bitflags0 | Particle_C.TRAJECTORY_FLAG
        cell_index = Mesh_C.NO_CELL
        unique_ID = Particle_C.UNIQUE_ID_COUNTER; Particle_C.UNIQUE_ID_COUNTER += 1
        crossings = 0

        # Trim the number of coordinates here to match
        # "position_coordinates" variable in ParticleInput_C
        p0 = (x0,y0, x0,y0, ux0,uy0, weight0, bitflags0, cell_index, unique_ID, crossings)

        # Second particle
        # An electron along the 22.5-degree line
        (x1, y1, z1) = (4.75, 0.0, 0.0)
        # Rotate through theta = pi/8 in x-y plane
# test
#        theta1 = math.pi/4.0
        theta1 = math.pi/8.0
        xnew = x1*np_m.cos(theta1)-y1*np_m.sin(theta1)
        ynew = x1*np_m.sin(theta1)+y1*np_m.cos(theta1)
        (x1, y1, z1) = (xnew, ynew, z1)
        (ux1, uy1, uz1) = (0.0, 0.0, 1000.0)
        weight1 = 2.0
        bitflags1 = 0b00 # bit flags variable is all zeroes
        # Turn trajectory flag ON.
        bitflags1 = bitflags1 | Particle_C.TRAJECTORY_FLAG
#        cell_index = -1
        cell_index = Mesh_C.NO_CELL
        unique_ID = Particle_C.UNIQUE_ID_COUNTER; Particle_C.UNIQUE_ID_COUNTER += 1
        crossings = 0
        
        # Trim the number of coordinates here to match
        # "position_coordinates" variable in ParticleInput_C
        p1 = (x1,y1, x1,y1, ux1,uy1, weight1, bitflags1, cell_index, unique_ID, crossings)

        particle_list = (p0, p1)
#        particle_list = (p1,)
        np = len(particle_list)

# Unit conversion                       
        return np, particle_list

    @staticmethod
    def test_electrons(type):
        """These electrons are specified as lists of coordinates and weights.
        """
        # Check that the caller knows this function returns 'listed' particles:
        if type != 'listed':
            error_msg = "test_electrons() called with type '" + str(type) + "' but it is of type 'listed'"
            sys.exit(error_msg)

        # Set initial phase-space coordinates and macroparticle weight for each particle

        # An electron at the y = 0 boundary

#        (x0, y0, z0) = (1.1, 0.1, 0.0)
#        (x0, y0, z0) = (1.0, 0.0, 0.0)
        (x0, y0, z0) = (4.75, 0.0, 0.0)
# Rotate through theta in x-y plane
        theta0 = math.pi/4.0
        xnew = x0*np_m.cos(theta0)-y0*np_m.sin(theta0)
        ynew = x0*np_m.sin(theta0)+y0*np_m.cos(theta0)
        (x0, y0, z0) = (xnew, ynew, z0)

#        print 'x0, y0, z0 =', x0, y0, z0

#        (ux0, uy0, uz0) = (3000.0, 2000.0, 1000.0)
        (ux0, uy0, uz0) = (0.0, 0.0, 1000.0)
        weight0 = 1.0 # number of electrons per macroparticle
        bitflags0 = 0b00 # bit flags variable is all zeroes
        # Turn on trajectory flag.
        bitflags0 = bitflags0 | Particle_C.TRAJECTORY_FLAG
        cell_index = Mesh_C.NO_CELL
        unique_ID = Particle_C.UNIQUE_ID_COUNTER; Particle_C.UNIQUE_ID_COUNTER += 1
        crossings = 0
        
        # Trim the number of coordinates here to match
        # "position_coordinates" variable in ParticleInput_C
        p0 = (x0,y0, x0,y0, ux0,uy0, weight0, bitflags0, cell_index, unique_ID, crossings)

        # Second particle
        (x1, y1, z1) = (1.1, 0.0, 10.0)
# Rotate through theta in x-y plane
        theta1 = math.pi/4.0
        xnew = x1*np_m.cos(theta1)-y1*np_m.sin(theta1)
        ynew = x1*np_m.sin(theta1)+y1*np_m.cos(theta1)
        (x1, y1, z1) = (xnew, ynew, z1)

#        print 'x1, y1, z1 =', x1, y1, z1

#        (ux1, uy1, uz1) = (3000.0, 2000.0, 1000.0)
        (ux1, uy1, uz1) = (0.0, 0.0, 1000.0)
        weight1 = 2.0
        bitflags1 = 0b00 # bit flags variable
#        TRAJECTORY_FLAG = 0b1
        # Turn on trajectory flag.
        bitflags1 = bitflags1 | Particle_C.TRAJECTORY_FLAG 
        unique_ID = Particle_C.UNIQUE_ID_COUNTER; Particle_C.UNIQUE_ID_COUNTER += 1
        crossings = 0
        
        # Trim the number of coordinates here to to match
        # "position_coordinates" variable in ParticleInput_C
        p1 = (x1,y1, x1,y1, ux1,uy1, weight1, bitflags1, cell_index, unique_ID, crossings)

        particle_list = (p0, p1)
#        particle_list = (p0,)
        np = len(particle_list)

# Unit conversion                       
#        velocity = [v*Convert.m_per_s for v in velocity]
        return np, particle_list
#     def test_electrons(type):ENDDEF

    @staticmethod
    def two_electrons(type):
        """These electrons are specified as lists of coordinates and weights.
        """
        # Check that the caller knows this function returns 'listed' particles:
        if type != 'listed':
            error_msg = "two_electrons() called with type '" + str(type) + "' but it is of type 'listed'"
            sys.exit(error_msg)

        # Set initial phase-space coordinates and macroparticle weight for each particle

        # An electron at the y = 0 boundary

        (x0, y0, z0) = (4.75, 0.0, 0.0)
# Rotate through theta in x-y plane
        theta0 = math.pi/4.0
        xnew = x0*np_m.cos(theta0)-y0*np_m.sin(theta0)
        ynew = x0*np_m.sin(theta0)+y0*np_m.cos(theta0)
        (x0, y0, z0) = (xnew, ynew, z0)

        print("First particle: x0, y0, z0 =", x0, y0, z0)

        (ux0, uy0, uz0) = (0.0, 0.0, 1000.0)
        weight0 = 1.0 # number of electrons per macroparticle
        bitflags0 = 0b00 # bit flags variable is all zeroes
        # Turn on trajectory flag.
        bitflags0 = bitflags0 | Particle_C.TRAJECTORY_FLAG
        cell_index = Mesh_C.NO_CELL
        unique_ID = Particle_C.UNIQUE_ID_COUNTER; Particle_C.UNIQUE_ID_COUNTER += 1
        crossings = 0
        
        # Trim the number of coordinates here to match
        # "position_coordinates" variable in ParticleInput_C
        p0 = (x0,y0, x0,y0, ux0,uy0, weight0, bitflags0, cell_index, unique_ID, crossings)

        # Second particle
        (x1, y1, z1) = (1.1, 0.0, 10.0)
# Rotate through theta in x-y plane
        theta1 = math.pi/4.0
        xnew = x1*np_m.cos(theta1)-y1*np_m.sin(theta1)
        ynew = x1*np_m.sin(theta1)+y1*np_m.cos(theta1)
        (x1, y1, z1) = (xnew, ynew, z1)

        print("Second particle: x1, y1, z1 =", x1, y1, z1)

        (ux1, uy1, uz1) = (0.0, 0.0, 1000.0)
        weight1 = 2.0
        bitflags1 = 0b00 # bit flags variable
#        TRAJECTORY_FLAG = 0b1
        # Turn on trajectory flag.
        bitflags1 = bitflags1 | Particle_C.TRAJECTORY_FLAG 
        unique_ID = Particle_C.UNIQUE_ID_COUNTER; Particle_C.UNIQUE_ID_COUNTER += 1
        crossings = 0
        
        # Trim the number of coordinates here to to match
        # "position_coordinates" variable in ParticleInput_C
        p1 = (x1,y1, x1,y1, ux1,uy1, weight1, bitflags1, cell_index, unique_ID, crossings)

        particle_list = (p0, p1)
#        particle_list = (p0,)
        np = len(particle_list)

# Unit conversion                       
#        velocity = [v*Convert.m_per_s for v in velocity]
        return np, particle_list
#     def two_electrons(type):ENDDEF

#class UserParticleDistributions_C(object):ENDCLASS

#class UserParticleMeshBoundaryConditions_C(ParticleBoundaryConditions_C):
#class UserParticleBoundaryFunctions_C(Particle_C):
#STARTCLASS
class UserParticleBoundaryFunctions_C(object):
    """UserParticleBoundaryFunctions_C implements callback functions
       (boundary conditions) for kinetic particles crossing marked
       mesh facets.

       See Particle_Module::ParticleMeshBoundaryConditions_C for naming
       scheme.
    """

    def __init__(self, position_coordinates, dx):
        """Initialize a UserParticleBoundaryFunctions_C object.

            :param position_coordinates: Example: ['x', 'y',]
            :param dx: The particle move vector
        
        """

        # Make aliases for quantities contained in particle_P that are needed to
        # implement various boundary conditions.
        self.position_coordinates = position_coordinates
        self.particle_dimension = len(self.position_coordinates)
        self.dx = dx

        # Create scratch
        precision = self.dx.dtype
        self.pcoord = np_m.empty(self.particle_dimension, dtype=precision)
        self.pvel = np_m.empty(self.particle_dimension, dtype=precision)

        return
#    def __init__(self, particle_P):ENDDEF

#class UserParticleBoundaryFunctions_C(object):
    @staticmethod
    def default_bc(p, species_name, facet_index, dx_fraction=None, facet_normal=None):
        """Default boundary condition for all species on all boundaries.
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print("Called", fncName)

        # Set the delete flag
        p['bitflags'] = p['bitflags'] | Particle_C.DELETE_FLAG

        return
#    def default_bc(p, species_name, facet_index, dx_fraction=None, facet_normal=None):ENDDEF


#class UserParticleBoundaryFunctions_C(object):
    @staticmethod
    def default_bc_at_rmin(p, species_name, facet_index, dx_fraction=None, facet_normal=None):
        """Default boundary condition for all particles incident on rmin.
        """

        printInfoInvoked = False
        
        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():'

        if printInfoInvoked is True:
            print("DnT INFO: %s Invoked by particle %s of species %s" % (fncName, p, species_name))

        # Set the delete flag
        p['bitflags'] = p['bitflags'] | Particle_C.DELETE_FLAG

        # Count the number/charge/energy of deleted particles

        return
#    def default_bc_at_rmin(p, species_name, facet_index, dx_fraction=None, facet_normal=None):ENDDEF
    
#class UserParticleBoundaryFunctions_C(object):
    def bc_at_rmin_for_test_electrons(self, p, species_name, facet_index, dx_fraction=None, facet_normal=None):
        """Boundary condition for test_electrons species incident on rmin.


           :param p: A full particle record
           :param str species_name: This is redundant since the function contains the
                                    name of the species, but may be useful for
                                    indexing.
           :param int facet_index: The mesh index of the facet that generated the call
                                  to this function.
           :param double dx_fraction: The fraction of the move vector traveled before
                                      the facet was crossed.
           :param double[] facet_normal: The unit vector normal to the facet crossed.

        """

        printInfoInvoked = False
        
        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():'

        if printInfoInvoked is True:
            print("DnT INFO: %s Invoked by particle %s of species %s" % (fncName, p, species_name))
        pDim = self.particle_dimension

        # Scratch space
        # self.pcoord can hold: x,y,z, (or subset)
        dxOutside = self.pcoord
        # self.pvel can hold: ux,uy,uz, (or subset)
        u = self.pvel

        # This is a reference to the actual move vector (see __init__() above)
        dx = self.dx # The the move-vector
       
        # Compute the move-vector past the reflecting surface
        for i in range(pDim):
            dxOutside[i] = (1.0 - dx_fraction)*dx[i]
        # Compute the component in the direction normal to the reflecting surface
        n_dot_dxOutside = np_m.dot(facet_normal, dxOutside)

        # Compute the component of velocity normal to the reflecting surface
        for i in range(pDim):
            u[i] = p[i+2*pDim]
        n_dot_u = np_m.dot(facet_normal, u)

        # Reflect the particle position in the surface to get its new final
        # position. Reflect the particle velocity too.
        i=0
        for coord in self.position_coordinates:
            p[coord] -= 2.0 * n_dot_dxOutside * facet_normal[i]
            ucomp = 'u'+coord
            p[ucomp] -= 2.0 * n_dot_u * facet_normal[i]
            i+=1
        
        # Count the number/charge/energy of deleted particles
        
        # Set the delete flag
#        p['bitflags'] = p['bitflags'] | Particle_C.DELETE_FLAG

        return
#    def bc_at_rmin_for_test_electrons(self, p, species_name, facet_index, dx_fraction=None, facet_normal=None):ENDDEF

#class UserParticleBoundaryFunctions_C(object):ENDCLASS
