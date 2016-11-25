# UserParticles Module

"""
These are static (like global) functions. No SELF variable!
"""
__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
__all__ = ['ParticleDistributions_C.drifting_Hplus',
           ]

import sys
import math
import numpy as np_M

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

class ParticleDistributions_C(object):
    """ParticleDistributions_C is to be edited by the user to
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

        # Trim the number of coordinates here if needed to match "position_coordinates" variable in DTparticleInput_C
        p0 = (x0,y0,z0, x0,y0,z0, ux0,uy0,uz0, weight0, bitflags0, cell_index)

        # Second particle
        # Particle moves diagonally across the mesh:
        (x1, y1, z1) = (9.5, 9.5, 9.5)
        (ux1, uy1, uz1) = (-2.0, -2.0, -2.0)
        weight1 = 3.0
        bitflags1 = 0b0 # bit flags variable
        # Turn ON trajectory flag.
        bitflags1 = bitflags1 | Particle_C.TRAJECTORY_FLAG 

        # Trim the number of coordinates here to match "position_coordinates" variable in DTparticleInput_C
        p1 = (x1,y1,z1, x1,y1,z1, ux1,uy1,uz1, weight1, bitflags1, cell_index)

#        particle_list = (p0,) # Need the final comma , for one particle.
#        particle_list = (p1,) # Need the final comma , for one particle.
        particle_list = (p0, p1)

        np = len(particle_list)

        return np, particle_list
#    def neutral_H(type):ENDDEF

#class ParticleDistributions_C(object):ENDCLASS
