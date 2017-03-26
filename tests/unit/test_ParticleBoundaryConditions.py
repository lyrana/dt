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

from UserUnits_Module import MyPlasmaUnits_C

class TestParticleBoundaryConditions(unittest.TestCase):
    """Test user-specified boundary conditions on particles.
    """
    
    def setUp(self):

        # Initialization code common to the tests go here...

        return

#class TestParticleBoundaryConditions(unittest.TestCase):
    def test_1_xy_absorbing_boundary(self):
        """ Check that particles are deleted correctly when they
            strike an absorbing boundary on a 2D Cartesian mesh.
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest: ', fncName

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
        pinCI.position_coordinates = ['x', 'y', 'z'] # determines the particle-storage dimensions
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
        pinCI.user_particles_module = 'UserParticles_3D'
        UPrt_M = im_M.import_module(pinCI.user_particles_module)
        pinCI.user_particles_class = UPrt_M.UserParticleDistributions_C

        # Make the particle object using pinCI input
        particleCI = Particle_C(pinCI, printFlag=False)

        ### Create a trajectory object and add it to particleCI

        # Create input object for trajectories
        trajinCI = DTtrajectoryInput_C()

        trajinCI.maxpoints = None # Set to None to get every point

        # Specify which particle variables to save.  This has the
        # form of a numpy dtype specification.
        trajinCI.explicitDict = {'names': ['x', 'ux', 'y', 'uy', 'Ex', 'Ey'], 'formats': [numpy.float32]*6}
        trajinCI.implicitDict = {'names': ['x', 'ux', 'phi'], 'formats': [numpy.float32]*3}
        trajinCI.neutralDict = {'names': ['x', 'ux', 'y', 'uy'], 'formats': [numpy.float32]*4}

        # Add a trajCI reference to the particle object
        pCI = particleCI # abbreviation
        pCI.trajCI = Trajectory_C(trajinCI, ctrlCI, pCI.explicit_species, pCI.implicit_species, pCI.neutral_species)

        # Initialize the particles
        printFlags = {}
        for sp in particleCI.species_names: printFlags[sp] = False
        particleCI.initialize_distributions(printFlags)

#  Mesh input for the particle mesh, including particle boundary conditions.

        # Here's the mesh definition for this test
        #from UserMesh_y_Fields_FE2D_Module import UserMesh_C
        from UserMesh_FE_XYZ_Module import UserMesh_C

        # Create a 2D Cartesian mesh to use for advancing the particles.  The particles
        # themselves are given 3D coordinates.

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
        # xmin_indx = '1'
        # xmax_indx = '2'
        # ymin_indx = '4'
        # ymax_indx = '8'
# This could be automated, given a list of the boundary names: ['xmin', 'xmax', ...]
        xminIndx = 1
        xmaxIndx = 2
        yminIndx = 4
        ymaxIndx = 8
        particleBoundaryDict = {'xmin': xminIndx,
                                'xmax': xmaxIndx,
                                'ymin': yminIndx,
                                'ymax': ymaxIndx,
                                }

        mi2DCI.particle_boundary_dict = particleBoundaryDict

        pmeshCI = UserMesh_C(mi2DCI, compute_dictionaries=True, compute_tree=True, plot_flag=False)
        # Add this to the particle object:
        particleCI.pmeshCI = pmeshCI

# Particle boundary-conditions

        # UserParticleMeshFunctions_C is where the facet-crossing callback
        # functions are defined.
#        user_particle_bcs_class = UPrt_M.UserParticleMeshBoundaryConditions_C
        userPMeshFnsClass = UPrt_M.UserParticleMeshFunctions_C # abbreviation

        # Make the particle-mesh boundary-conditions object and add it
        # to the particle object.  The user has to supply the
        # facet-crossing callback functions in the
        # UserParticleMeshFunctions_C object above.

        spNames = particleCI.species_names
        pmeshBCCI = ParticleMeshBoundaryConditions_C(spNames, pmeshCI, userPMeshFnsClass, printFlag=False)
        particleCI.pmesh_bcCI = pmeshBCCI

# Get the initial cell index of each particle.

# Should this be something the pmesh computes?  No: pmesh computes the
# index of a single particle.  It doesn't know the particle storage
# infrastructure.
        particleCI.compute_mesh_cell_indices()

# Advance the particles for nsteps

        print "Moving", pCI.get_total_particle_count(), "particles for", ctrlCI.nsteps, "timesteps"
        for istep in xrange(ctrlCI.nsteps):
            particleCI.move_neutral_particles(ctrlCI.dt)

        return
#    def test_2D_xy_absorbing_boundary(self):ENDDEF

#class TestParticleBoundaryConditions(unittest.TestCase):ENDCLASS

#class TestParticleBoundaryConditions(unittest.TestCase):
    def test_2_r_theta_absorbing_boundary(self):
        """ Check that particles are deleted correctly when they
            strike an absorbing boundary.
        """

        testName = sys._getframe().f_code.co_name
        fncName = '('+__file__+') ' + testName + '():\n'
        print '\ntest: ', fncName

        # if os.environ.get('DISPLAY') is None:
        #     plotFlag=False
        # else:
        #     plotFlag=True

        ### Particle species input

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

        # Make the particle object using pinCI input
        particleCI = Particle_C(pinCI, printFlag=False)

        ### Add a ref to a Trajectory_C object to particleCI

        # Create input object for trajectories
        trajinCI = DTtrajectoryInput_C()

        trajinCI.maxpoints = None # Set to None to get every point

        # Specify which particle variables to save.  This has the
        # form of a numpy dtype specification.
        trajinCI.explicitDict = {'names': ['x', 'ux', 'y', 'uy', 'Ex', 'Ey'], 'formats': [numpy.float32]*6}
        trajinCI.implicitDict = {'names': ['x', 'ux', 'phi'], 'formats': [numpy.float32]*3}
        trajinCI.neutralDict = {'names': ['x', 'ux', 'y', 'uy'], 'formats': [numpy.float32]*4}

        # # Initialize the particles
        # printFlags = {}
        # for sp in particleCI.species_names: printFlags[sp] = False
        # particleCI.initialize_distributions(printFlags)

        ###  Mesh and Fields input for the particle mesh.

        # Here's the mesh definition for this test
        from UserMesh_y_Fields_FE2D_Module import UserMesh_C

        ## Get the mesh from an existing file

        miCI = DTmeshInput_C()
        miCI.mesh_file = 'quarter_circle_mesh_crossed.xml'
        miCI.particle_boundary_file='Pbcs_quarter_circle_mesh_crossed.xml'
        # These are the boundary-name -> int pairs used to mark mesh facets:
        rminIndx = 1
        rmaxIndx = 2
        thminIndx = 4
        thmaxIndx = 8
        particleBoundaryDict = {'rmin': rminIndx,
                                'rmax': rmaxIndx,
                                'thmin': thminIndx,
                                'thmax': thmaxIndx,
                                }

        miCI.particle_boundary_dict = particleBoundaryDict

        pmeshCI = UserMesh_C(miCI, compute_dictionaries=True, compute_tree=True, plot_flag=False)
        # Add this to the particle object:
        particleCI.pmeshCI = pmeshCI

        ## Get the electric field from an existing file

        # The following value should correspond to the element degree
        # used in the potential from which negE was obtained
        phiElementDegree = 1

        if phiElementDegree == 1:
            # For linear elements, grad(phi) is discontinuous across
            # elements. To represent this field, we need Discontinuous Galerkin
            # elements.
            electricFieldElementType = 'DG'
        else:
            electricFieldElementType = 'Lagrange'

        # Create the negative electric field directly on the particle mesh
        negElectricField = Field_C(meshCI=particleCI.pmeshCI,
                                          element_type=electricFieldElementType,
                                          element_degree=phiElementDegree-1,
                                          field_type='vector')

        file = df_M.File('negE2D_crossed.xml')
        file >> negElectricField.function

        ## Particle boundary-conditions

        # UserParticleMeshFunctions_C is where the facet-crossing callback
        # functions are defined.
        userPMeshFnsClass = UPrt_M.UserParticleMeshFunctions_C # abbreviation

        # Make the particle-mesh boundary-conditions object and add it
        # to the particle object.
        spNames = particleCI.species_names
        pmeshBCCI = ParticleMeshBoundaryConditions_C(spNames, pmeshCI, userPMeshFnsClass, printFlag=False)
        particleCI.pmesh_bcCI = pmeshBCCI


        ctrlCI = DTcontrol_C()

        # These are fast electrons, so the timestep is small
        ctrlCI.dt = 1.0e-6
        ctrlCI.nsteps = 13

        # The trajectory object can now be created and added to particleCI
        pCI = particleCI
        pCI.trajCI = Trajectory_C(trajinCI, ctrlCI, pCI.explicit_species, pCI.implicit_species, pCI.neutral_species)


        # Initialize the particles
        printFlags = {}
        for sp in pCI.species_names: printFlags[sp] = False
        pCI.initialize_distributions(printFlags)

        # Get the initial cell index of each particle.
        pCI.compute_mesh_cell_indices()

        ### Particle loop

        print "Moving", pCI.get_total_particle_count(), "particles for", ctrlCI.nsteps, "timesteps"

        for istep in xrange(ctrlCI.nsteps):

            if pCI.trajCI is not None:
#                print 'pCI.trajCI.skip:', pCI.trajCI.skip
                if istep % pCI.trajCI.skip == 0:
                    pCI.record_trajectory_data(neg_E_field=negElectricField)

            pCI.move_particles_in_electrostatic_field(ctrlCI.dt, negElectricField)

        # Record the LAST point on the particle trajectory
        if pCI.trajCI is not None:
                pCI.record_trajectory_data(neg_E_field=negElectricField)

        # Plot the trajectory onto the particle mesh

        mesh = pCI.pmeshCI.mesh
        holdPlot = True # Set to True to stop the plot from disappearing.
        pCI.trajCI.plot_trajectories_on_mesh(mesh, testName, hold_plot=holdPlot) # Plots trajectory spatial coordinates on top of the particle mesh

        return
#    def test_2D_r_theta_absorbing_boundary(self):

#class TestParticleBoundaryConditions(unittest.TestCase):ENDCLASS

if __name__ == '__main__':
    unittest.main()
