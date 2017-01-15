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

from UserMesh_y_Fields_FE2D_Module import UserMesh_C

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

        self.ctrlCI = DTcontrol_C()

        self.ctrlCI.dt = 1.0e-6
        self.ctrlCI.nsteps = 13

#
# Particle species input
#

        # Create an instance of the DTparticleInput class
        pinCI = DTparticleInput_C()
        # Initialize particles
        pinCI.precision = numpy.float64
        pinCI.particle_integration_loop = 'loop-on-particles'
        pinCI.position_coordinates = ['x', 'y',] # determines the particle-storage dimensions
        pinCI.force_components = ['x', 'y',]
        pinCI.force_precision = numpy.float64

        # Specify the particle-species properties
        # 1. electrons
        pinCI.particle_species = (('trajelectrons',
                             {'initial_distribution_type' : 'listed',
                              'charge' : -1.0*MyPlasmaUnits_C.elem_charge,
                              'mass' : 1.0*MyPlasmaUnits_C.electron_mass,
                              'dynamics' : 'explicit',
#                              'dynamics' : 'implicit',
                              }
                             ),
                            )

        # Provide the particle distributions for the above species
        pinCI.user_particles_module = 'UserParticles_2D_e'
        UPrt_M = im_M.import_module(pinCI.user_particles_module)
        pinCI.user_particles_class = UPrt_M.UserParticleDistributions_C

        # Particle boundary-conditions

# Add a ref to a ParticleBoundaryConditions_C object to particleCI
        pinCI.user_particle_bcs_class = UPrt_M.UserParticleBoundaryConditions_C
#
#  Mesh and Fields input for the particle mesh.
#
        # Make the mesh, boundary-conditions markers & fields from saved files
# 'crossed' diagonals
#        self.pmesh2DCI = Mesh_C(meshFile='quarter_circle_mesh_crossed.xml', computeDictionaries=True, computeTree=True, plotFlag=False)
# Replace with UserMesh_C because this includes a boundary_marker.
#        self.pmesh2DCI = UserMesh_C(meshFile='quarter_circle_mesh_crossed.xml', computeDictionaries=True, computeTree=True, plotFlag=False)

# 'left' diagonal
#        self.pmesh2DCI = Mesh_C(meshFile='quarter_circle_mesh_left.xml', computeDictionaries=True, computeTree=True, plotFlag=False)

        # Get the mesh from an existing file

        miCI = DTmeshInput_C()
        miCI.mesh_file = 'quarter_circle_mesh_crossed.xml'
        miCI.particle_boundary_file='Pbcs_quarter_circle_mesh_crossed.xml'
        # These are the boundary-name -> int pairs used to mark mesh facets:
        rmin_indx = 1
        rmax_indx = 2
        thmin_indx = 4
        thmax_indx = 8
        particleBoundaryDict = {'rmin': rmin_indx,
                                'rmax': rmax_indx,
                                'thmin': thmin_indx,
                                'thmax': thmax_indx,
                                }

        miCI.particle_boundary_dict = particleBoundaryDict


#        particleCI.pmeshCI.particle_boundary_dict = particleBoundaryDict

        # Add particle mesh to pinCI
#        pinCI.pmeshCI = UserMesh_C(meshFile='quarter_circle_mesh_crossed.xml', particleBoundaryFile='Pbcs_quarter_circle_mesh_crossed.xml', computeDictionaries=True, computeTree=True, plotFlag=False)
        pinCI.pmeshCI = UserMesh_C(miCI, computeDictionaries=True, computeTree=True, plotFlag=False)

        # Make the particle object from pinCI
        particleCI = Particle_C(pinCI, printFlag=False)

        # The following value should correspond to the element degree
        # used in the potential from which negE was obtained
        phi_element_degree = 1

        if phi_element_degree == 1:
            # For linear elements, grad(phi) is discontinuous across
            # elements. To represent this field, we need Discontinuous Galerkin
            # elements.
            electric_field_element_type = 'DG'
        else:
            electric_field_element_type = 'Lagrange'

        # Create the negative electric field directly on the particle mesh
        self.neg_electric_field = Field_C(meshCI=particleCI.pmeshCI,
                                          element_type=electric_field_element_type,
                                          element_degree=phi_element_degree-1,
                                          field_type='vector')

        file = df_M.File('negE2D_crossed.xml')
        file >> self.neg_electric_field.function

        # Mark these subdomains (boundaries) with non-zero integers
        # Name the boundaries used to apply particle
        # boundary-conditions. Integers will be assigned to them in
        # UserMesh_C.
#        particleBoundaryNames = ['rmin', 'rmax', 'thmin', 'thmax']
#        pinCI.particle_boundary_names = ['rmin', 'rmax', 'thmin', 'thmax']

        # Make the particle boundary-conditions object and add it to the particle object

        pmbcCI = ParticleMeshBoundaryConditions_C(pinCI, particleCI, printFlag=False)
        particleCI.pbcCI = pmbcCI

# Add a ref to a Trajectory_C object to particleCI

        # Create input object for trajectories
        self.trajinCI = DTtrajectoryInput_C()

        self.trajinCI.maxpoints = None # Set to None to get every point

        # Specify which particle variables to save.  This has the
        # form of a numpy dtype specification.
        self.trajinCI.explicitDict = {'names': ['x', 'ux', 'y', 'uy', 'Ex', 'Ey'], 'formats': [numpy.float32]*6}
        self.trajinCI.implicitDict = {'names': ['x', 'ux', 'phi'], 'formats': [numpy.float32]*3}
        self.trajinCI.neutralDict = {'names': ['x', 'ux', 'y', 'uy'], 'formats': [numpy.float32]*4}

        # Add reference to particleCI
        pCI = particleCI
        pCI.trajCI = Trajectory_C(self.trajinCI, self.ctrlCI, pCI.explicit_species, pCI.implicit_species, pCI.neutral_species)


#    def test_2D_absorbing_boundary(self):ENDDEF

#class TestParticleBoundaryConditions(unittest.TestCase):ENDCLASS

if __name__ == '__main__':
    unittest.main()
