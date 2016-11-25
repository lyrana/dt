# UserParticles Module

"""
These are static (like global) functions. No SELF variable!
"""
__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
__all__ = ['ParticleDistributions_C.initial_electrons',
           'ParticleDistributions_C.initial_ions',
           'ParticleDistributions_C.few_initial_electrons',
           'ParticleDistributions_C.few_initial_ions',
           ]

import sys
from Particles_Module import Particles_C
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

#ClassClassClassClassClassClassClassClassclass
class ParticleDistributions_C(object):
    """ParticleDistributions_C is to be edited by the user to
       specify initial particle distributions.  The units are MKS by
       default (i.e., if no conversion factor is applied), but any
       units can be used provided the conversion to MKS is available
       in the UserUnits_M.py module.
    """

    # The spatial coordinates of all particles
#    position_coordinates = ['x', 'y',]


# Initialize a few particles from each species
# 
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
        (x0, y0, z0) = (1.0, 0.0, 0.0)
        (ux0, uy0, uz0) = (3000.0, 2000.0, 1000.0)
        weight0 = 1.0 # number of electrons per macroparticle
        bitflags0 = 0b00 # bit flags variable
        bitflags0 = bitflags0 | Particles_C.TRAJECTORY_FLAG # use low bit for trajectory flag.

        # Trim the number of coordinates here to match "position_coordinates" class variable above
        p0 = (x0,y0, ux0,uy0, weight0, bitflags0)

        # Image particle with different weight
        (x1, y1, z1) = (4.75, 0.0, 0.0)
        (ux1, uy1, uz1) = (3000.0, 2000.0, 1000.0)
        weight1 = 2.0
        bitflags1 = 0b00 # bit flags variable
#        TRAJECTORY_FLAG = 0b1
        bitflags1 = bitflags1 | Particles_C.TRAJECTORY_FLAG # set trajectory bit to 1.

        # Trim the number of coordinates here to match "position_coordinates" class variable above
        p1 = (x1,y1, ux1,uy1, weight1, bitflags1)

        particle_list = (p0, p1)
        np = len(particle_list)

# Unit conversion                       
#        velocity = [v*Convert.m_per_s for v in velocity]
        return np, particle_list

