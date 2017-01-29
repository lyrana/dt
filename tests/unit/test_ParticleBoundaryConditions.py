#!/usr/bin/env python

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
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

# Here's the mesh definition for this test
#from UserMesh_y_Fields_FE2D_Module import UserMesh_C
from UserMesh_FE_XYZ_Module import UserMesh_C

from UserUnits_Module import MyPlasmaUnits_C

class TestParticleBoundaryConditions(unittest.TestCase):
    """Test user-specified boundary conditions on particles.
    """
    
    def setUp(self):

        # Initialization code common to the tests go here...

        return

#class TestParticleBoundaryConditions(unittest.TestCase):
    def test_2D_absorbing_boundary(self):
        """ Check that particles are deleted correctly when they
            strike an absorbing boundary.
        """

        fncname = sys._getframe().f_code.co_name
        print '\ntest: ', fncname, '('+__file__+')'

        # if os.environ.get('DISPLAY') is None:
        #     plotFlag=False
        # else:
        #     plotFlag=True

        ctrlCI = DTcontrol_C()

        ctrlCI.dt = 0.5
        ctrlCI.nsteps = 100

# Particle species input

        # Create an instance of the DTparticleInput class
        pinCI = DTparticleInput_C()
        # Initialize particles
        pinCI.precision = numpy.float64
        pinCI.particle_integration_loop = 'loop-on-particles'
 #       pinCI.position_coordinates = ['x', 'y',] # determines the particle-storage dimensions
        pinCI.position_coordinates = ['x', 'y', 'z'] # determines the particle-storage dimensions
#        pinCI.force_components = ['x', 'y',]
        pinCI.force_precision = numpy.float64

        # Specify the particle-species properties
        # 1. electrons
        pinCI.particle_species = (('neutral_H',
                             {'initial_distribution_type' : 'listed',
                              'charge' : 0.0,
                              'mass' : 1.0*MyPlasmaUnits_C.AMU,
                              'dynamics' : 'explicit',
#                              'number_per_cell' : 12,
                              }
                             ),
                            )

        # Provide the particle distributions for the above species
        # This could be done more like how the mesh is specified:
        # import UserParticles_3D as UPrt_M
#        pinCI.user_particles_module = 'UserParticles_2D_e'
        pinCI.user_particles_module = 'UserParticles_3D'
        UPrt_M = im_M.import_module(pinCI.user_particles_module)
        pinCI.user_particles_class = UPrt_M.UserParticleDistributions_C

        # Particle boundary-conditions

        # Add a ref to a ParticleBoundaryConditions_C object to
        # particleCI.  This is where the facet-crossing callback
        # functions are defined.
        pinCI.user_particle_bcs_class = UPrt_M.UserParticleMeshBoundaryConditions_C

#  Mesh input for the particle mesh.

        # Create a mesh to push the particles on.

        mi2DCI = DTmeshInput_C()

        # 2D mesh input
        mi2DCI = DTmeshInput_C()
        (xmin, ymin) = (-10.0, -10.0)
        (xmax, ymax) = ( 10.0,  10.0)
        mi2DCI.pmin = df_M.Point(xmin, ymin)
        mi2DCI.pmax = df_M.Point(xmax, ymax)
        mi2DCI.cells_on_side = (4, 2)
        # These are the (int boundary-name) pairs used to mark mesh
        # facets. The string value of the int is used as the index.
        xmin_indx = '1'
        xmax_indx = '2'
        ymin_indx = '4'
        ymax_indx = '8'
        # particleBoundaryDict = {'xmin': xmin_indx,
        #                         'xmax': xmax_indx,
        #                         'ymin': ymin_indx,
        #                         'ymax': ymax_indx,
        #                         }
        # Note the order, which is reversed from the fieldBoundaryDict
        # order.
        particleBoundaryDict = {xmin_indx:'xmin',
                                xmax_indx:'xmax',
                                ymin_indx:'ymin',
                                ymax_indx:'ymax',
                                }

        mi2DCI.particle_boundary_dict = particleBoundaryDict

        # Make the mesh and point pinCI.pmeshCI to it.
        pinCI.pmeshCI = UserMesh_C(mi2DCI, computeDictionaries=True, computeTree=True, plotFlag=False)

        # Make the particle object using pinCI input
        # (The pmeshCI ref gets copied from pinCI to particleCI)
        particleCI = Particle_C(pinCI, printFlag=False)

        # Initialize the particles
        printFlags = {}
        for sp in particleCI.species_names: printFlags[sp] = False
        particleCI.initialize_distributions(printFlags)
        # Get the initial cell index of each particle.
        particleCI.compute_mesh_cell_indices()


        # Make the particle-mesh boundary-conditions object and add it
        # to the particle object.  The user has to supply the
        # facet-crossing callback functions in the
        # UserParticleMeshBoundaryConditions_C object above.

        pmbcCI = ParticleMeshBoundaryConditions_C(pinCI, particleCI, printFlag=False)
        particleCI.pbcCI = pmbcCI

# Add a ref to a Trajectory_C object to particleCI

        # Create input object for trajectories
        trajinCI = DTtrajectoryInput_C()

        trajinCI.maxpoints = None # Set to None to get every point

        # Specify which particle variables to save.  This has the
        # form of a numpy dtype specification.
        trajinCI.explicitDict = {'names': ['x', 'ux', 'y', 'uy', 'Ex', 'Ey'], 'formats': [numpy.float32]*6}
        trajinCI.implicitDict = {'names': ['x', 'ux', 'phi'], 'formats': [numpy.float32]*3}
        trajinCI.neutralDict = {'names': ['x', 'ux', 'y', 'uy'], 'formats': [numpy.float32]*4}

        # Add reference to particleCI
        pCI = particleCI
        pCI.trajCI = Trajectory_C(trajinCI, ctrlCI, pCI.explicit_species, pCI.implicit_species, pCI.neutral_species)


        # Now advance the particles for nsteps
        print "Moving", pCI.get_total_particle_count(), "particles for", ctrlCI.nsteps, "timesteps"
        for istep in xrange(ctrlCI.nsteps):
            particleCI.move_neutral_particles(ctrlCI.dt)


#    def test_2D_absorbing_boundary(self):ENDDEF

#class TestParticleBoundaryConditions(unittest.TestCase):ENDCLASS

if __name__ == '__main__':
    unittest.main()
