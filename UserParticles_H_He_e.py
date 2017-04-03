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
from Particle_Module import Particle_C
from Dolfin_Module import Mesh_C
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

    # The spatial coordinates of all particles
# Moved this to ParticleInput_C:
#    position_coordinates = ['x', 'y', 'z']

# These distributions are constant-density in space
    @staticmethod
    def initial_electrons(p):
        density = 2.0e10*Convert.number_per_m3
        temperature = 2.0*Convert.eV
        velocity = [0.0, 10.0, 0.0]
        velocity = [v*Convert.m_per_s for v in velocity]
        return density, velocity, temperature

    @staticmethod
    def initial_ions(p):
        density = 2.0e10*Convert.number_per_m3
        temperature = 0.1*Convert.eV
        velocity = [0.0, 0.1, 0.0]
        velocity = [v*Convert.m_per_s for v in velocity]
        return density, velocity, temperature

# Initialize a few particles from each species
    @staticmethod
    def plasmaelectrons(type):
        """The plasma electrons are specified as lists of coordinates and weights.
        """
        # Check that the caller knows this function returns 'listed' particles:
        if type != 'listed':
            error_msg = "plasmaelectrons() called with type '" + str(type) + "' but it is of type 'listed'"
            sys.exit(error_msg)

        # Set initial phase-space coordinates and macroparticle weight for each particle

        (x0, y0, z0) = (0.0, 0.01, 0.02)
# test cell_index:
#        (x0, y0, z0) = (0.029, 0.016, 0.02)
        (ux0, uy0, uz0)  = (3000.0, 2000.0, 1000.0)

        weight0 = 1.0 # number of electrons per macroparticle
        bitflags0 = 0b00 # initialize all bits to 0
        bitflags0 = bitflags0 | Particle_C.TRAJECTORY_FLAG # turn on trajectory flag
#        cell_index = -1
        cell_index = Mesh_C.NO_CELL
#        print "plasmaele... Particle_C.TRAJECTORY_FLAG, bitflags0 = ", Particle_C.TRAJECTORY_FLAG, bitflags0

        p0 = (x0,y0,z0, x0,y0,z0, ux0,uy0,uz0, weight0,bitflags0,cell_index)

        # Image particle with different weight
        x1 = y1 = 0.0; z1 = -z0
#test cell_index
#        (x1, y1, z1) = (0.001, 0.029, 0.02)

        ux1 = uy1 = 0.0; uz1 = -uz0
        weight1 = 2.0
        bitflags1 = 0b00 # initialize all bits to 0
#        bitflags1 = bitflags0 | Particle_C.TRAJECTORY_FLAG # use low bit for trajectory flag.

        p1 = (x1,y1,z1, x1,y1,z1, ux1,uy1,uz1, weight1,bitflags1,cell_index)

        particle_list = (p0, p1)
        np = len(particle_list)

# Unit conversion                       
#        velocity = [v*Convert.m_per_s for v in velocity]
        return np, particle_list

# 

    @staticmethod
    def Hplus(type):
        """The Hplus ions are specified as lists of coordinates and weights.
        """
        # Check that the caller knows this function returns 'listed' particles:
        if type != 'listed':
            error_msg = "Hplus() called with type '" + str(type) + "' but it is of type 'listed'"
            sys.exit(error_msg)

        # Set phase-space coordinates and macroparticle weight for each particle

        (x0, y0, z0) = (0.0, 1.e-5, 2e-5)
# test cell_index
#        (x0, y0, z0) = (-0.001, 0.0, 2e-5)

#        ux0 = 0.3; uy0 = 0.2; uz0 = 0.1
        (ux0, uy0, uz0) = (0.3, 0.2, 0.1)
        weight0 = 3.0
        bitflags0 = 0b00 # initialize all bits to 0
        bitflags0 = bitflags0 | Particle_C.TRAJECTORY_FLAG # turn on trajectory flag
        cell_index = Mesh_C.NO_CELL
#        cell_index = -1

        # Trim the number of coordinates here to match "position_coordinates" variable in ParticleInput_C
        p0 = (x0,y0,z0, x0,y0,z0, ux0,uy0,uz0, weight0, bitflags0, cell_index)

        # Image particle with different weight
        x1 = -x0; y1 = -y0; z1 = -z0
        ux1 = -ux0; uy1 = -uy0; uz1 = -uz0
        weight1 = 4.0
        bitflags1 = 0b01 # initialize all bits to 0
#        bitflags1 = bitflags1 | Particle_C.TRAJECTORY_FLAG # use low bit for trajectory flag.

        # Trim the number of coordinates here to match "position_coordinates" variable in ParticleInput_C
        p1 = (x1,y1,z1, x1,y1,z1, ux1,uy1,uz1, weight1, bitflags1, cell_index)

        particle_list = (p0, p1)
        np = len(particle_list)

# Unit conversion                       
#        velocity = [v*Convert.m_per_s for v in velocity]
        return np, particle_list

    @staticmethod
    def one_electron(type):
        """Initialize one electron
        """
        # Check that the caller knows this function returns 'listed' particles:
        if type != 'listed':
            error_msg = "one_electron() called with type '" + str(type) + "' but it is of type 'listed'"
            sys.exit(error_msg)

        # Set initial phase-space coordinates and macroparticle weight for each particle

        (x, y, z) = (0.0, 0.01, 0.02)
        (ux, uy, uz)  = (3000.0, 2000.0, 1000.0)

        weight = 1.0 # number of electrons per macroparticle
        bitflags = 0b00 # initialize all bits to 0
        bitflags = bitflags | Particle_C.TRAJECTORY_FLAG # turn on trajectory flag
        cell_index = Mesh_C.NO_CELL
#        cell_index = -1

        p = (x,y,z, x,y,z, ux,uy,uz, weight,bitflags,cell_index)

        particle_list = (p,)
        np = len(particle_list)

# Unit conversion                       
#        velocity = [v*Convert.m_per_s for v in velocity]
        return np, particle_list

    def few_initial_electrons(p):
#        charge = 1.0*Convert.coulomb
        position = [1.0, 1.0, 1.0]
        velocity = [0.0, 10.0, 0.0]
        velocity = [v*Convert.m_per_s for v in velocity]
        return density, velocity, temperature

    @staticmethod
    def few_initial_ions(p):
        density = 2.0e10*Convert.number_per_m3
        temperature = 0.1*Convert.eV
        velocity = [0.0, 0.1, 0.0]
        velocity = [v*Convert.m_per_s for v in velocity]
# Convert velocity to p/m?
        return density, velocity, temperature

#class UserParticleDistributions_C(object):ENDCLASS

#from Particle_Module import ParticleBoundaryConditions_C as PBC_C

#class UserParticleBoundaryConditions_C(ParticleBoundaryConditions_C):
class UserParticleBoundaryConditions_C(object):
    """UserParticleBoundaryConditions_C is to be edited by the user to
       provide less commonly-used treatments of particles incident on a mesh boundary.
    """

    # Static class variables

    # Particle boundary-conditions are invoked by turning on bits.
    # Continue from the number of standard boundary-conditions.
#    bitOffset = ParticleBoundaryConditions_C.NUMBER_OF_STANDARD_BCs # Is the classname needed?
#    ISEE = 0b1 << bitOffset        # Ion-stimulated electron emission
#    SEE  = 0b1 << 1 + bitOffset    # Secondary-electron emission

    # The spatial coordinates of all particles
# Moved this to ParticleInput_C:
#    position_coordinates = ['x', 'y', 'z']

# These distributions are constant-density in space
#    @staticmethod
#    def initial_electrons(p):
#        density = 2.0e10*Convert.number_per_m3

#class UserParticleBoundaryConditions_C(object):ENDCLASS
