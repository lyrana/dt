#!/usr/bin/env python

__version__ = 0.1
__author__ = 'Copyright (C) 2017 L. D. Hughes'
#__all__ = []

import sys
import os
import numpy
import importlib as im_M
import unittest

import dolfin as df_M

from DT_Module import DTmeshInput_C
from DT_Module import DTparticleInput_C
from DT_Module import DTcontrol_C
from DT_Module import DTtrajectoryInput_C

#from Dolfin_Module import Mesh_C
from Dolfin_Module import Field_C

from SegmentedArrayPair_Module import SegmentedArray_C
from Particle_Module import Particle_C
from Particle_Module import ParticleMeshBoundaryConditions_C
from Trajectory_Module import Trajectory_C

from UserUnits_Module import MyPlasmaUnits_C

class TestParticleGeneration(unittest.TestCase):
    """Test user-specified boundary conditions on particles.
    """
    
    def setUp(self):

        # Initialization code common to the tests go here...

        # Common particle inputs
        self.pinCI = DTparticleInput_C()

        self.pinCI.precision = numpy.float64
        self.pinCI.particle_integration_loop = 'loop-on-particles'
        self.pinCI.force_precision = numpy.float64

        # Common mesh inputs
        self.miCI = DTmeshInput_C()

        return

#class TestParticleGeneration:
    def test_1D_particle_source_region(self):
        """Set up a one-dimensional source region and display it.

        """

        pinCI = self.pinCI

        pinCI.position_coordinates = ['x',] # determines the particle-storage dimensions
        pinCI.force_components = ['x',]

        ### Particle input

        # Specify the particle properties
        # 1. electrons
        pinCI.particle_species = (('testelectrons',
                             {'initial_distribution_type' : 'listed',
                              'charge' : -1.0*MyPlasmaUnits_C.elem_charge,
                              'mass' : 1.0*MyPlasmaUnits_C.electron_mass,
                              'dynamics' : 'explicit',
                              }
                             ),
                            )

        # Provide the particle distributions (lists of particles, functions, etc.)
        pinCI.user_particles_module ="UserParticles_1D"
        
#        UPrt_M = im_M.import_module(pinCI.user_particles_module)
#        pinCI.user_particles_class = UPrt_C = UPrt_M.UserParticleDistributions_C

        # Make the particle storage array for all species
#        particleCI = Particle_C(pinCI, printFlag=False)

        ### 1D mesh input

        mi1DCI = self.miCI

        mi1DCI = DTmeshInput_C()
        mi1DCI.pmin = df_M.Point(-10.0)
        mi1DCI.pmax = df_M.Point(10.0)
        mi1DCI.cells_on_side = (4)

        ## Input for particle-mesh boundary-conditions

        # These are the (int boundary-name) pairs used to mark mesh
        # facets.
        xmin_indx = 1
        xmax_indx = 2
        particleBoundaryDict = {'xmin': xmin_indx,
                                'xmax': xmax_indx,
                                }

        # Add these to the mesh input
        mi1DCI.particle_boundary_dict = particleBoundaryDict

        ## Input for particle-generation regions

        # These are the (int source-name) pairs used to mark mesh
        # cells. The string value of the int is used as the index.
        sourceX1_indx = '1'
        sourceX2_indx = '2'
        # Note the order, which is reversed from the fieldBoundaryDict
        # order.
        particleSourceDict = {sourceX1_indx:'sourceX1',
                                sourceX2_indx:'sourceX2',
                                }

        # Add these to the mesh input
        mi1DCI.particle_source_dict = particleSourceDict


        # Create the 1D particle mesh
#        self.pmesh1DCI = UserMesh_C(mi1DCI, computeDictionaries=True, computeTree=True, plotFlag=plotFlag)




        return

    def test_2D_particle_source_region(self):
        """Set up a two-dimensional source region and display it.

        """

        mi2DCI = self.miCI

        mi2DCI.pmin = df_M.Point(-10.0, -10.0)
        mi2DCI.pmax = df_M.Point(10.0, 10.0)
        mi2DCI.cells_on_side = (4, 2)

        ## Boundary conditions for the particles on this mesh

        # Create a 2D particle mesh
#        self.pmesh2DCI = UserMesh_C(mi2DCI, computeDictionaries=True, computeTree=True, plotFlag=plotFlag)

        return

#class TestParticleGeneration(unittest.TestCase):ENDCLASS

if __name__ == '__main__':
    unittest.main()
