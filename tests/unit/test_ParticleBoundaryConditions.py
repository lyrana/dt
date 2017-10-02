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

from DT_Module import DTcontrol_C

from Dolfin_Module import Field_C

from SegmentedArrayPair_Module import SegmentedArray_C

from Particle_Module import *
from Particle_Module import *

from Trajectory_Module import *

from UserUnits_Module import MyPlasmaUnits_C

#STARTCLASS
class TestParticleBoundaryConditions(unittest.TestCase):
    """Test user-specified boundary conditions on particles.
    """
    
    def setUp(self):

        # Initialization code common to the tests go here...

        return

#class TestParticleBoundaryConditions(unittest.TestCase):
    def test_1_2D_x_y_absorbing_boundary(self):
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
        ctrlCI.n_timesteps = 100

        # Create an instance of the DTparticleInput class
        pinCI = ParticleInput_C()
        # Settings common to all species
        pinCI.precision = numpy.float64
        pinCI.particle_integration_loop = 'loop-on-particles'
        pinCI.position_coordinates = ['x', 'y', 'z'] # determines the particle-storage dimensions
        pinCI.force_precision = numpy.float64

        # Specify the particle species properties

        # Specify the particle-species
        # 1. electrons
#         pinCI.particle_species = (('neutral_H',
#                              {'initial_distribution_type' : 'listed',
#                               'charge' : 0.0,
#                               'mass' : 1.0*MyPlasmaUnits_C.AMU,
#                               'dynamics' : 'explicit',
# #                              'number_per_cell' : 12,
#                               }
#                              ),
#                             )

        speciesName = 'neutral_H' # Use this name later to initialize this
                                  # species or to create a source of this
                                  # species.
        charge = 0.0
        mass = 1.0*MyPlasmaUnits_C.AMU
        dynamics = 'explicit'
        neutralHCI = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        # Add the neutral hydrogen to particle input
        pinCI.particle_species = (neutralHCI,
                                 )

        ## Make the particle storage array for all species.
        particleCI = Particle_C(pinCI, print_flag=False)

        # Provide the particle distributions for the above species
        # This could be done more like how the mesh is specified:
        # import UserParticles_3D as UPrt_M
        # Particles have a 3D/3D phase-space even though mesh is just 2D
        userParticleModule = 'UserParticles_3D'

        # Import this module
        UPrt_M = im_M.import_module(userParticleModule)

        particleCI.user_particle_module = userParticleModule
        particleCI.user_particle_class = userParticleClass = UPrt_M.UserParticleDistributions_C

        ### Create a trajectory object and add it to particleCI

        # Create input object for trajectories
        trajinCI = TrajectoryInput_C()

        trajinCI.maxpoints = None # Set to None to get every point

        # Specify which particle variables to save.  This has the
        # form of a numpy dtype specification.
        trajinCI.explicit_dict = {'names': ['x', 'ux', 'y', 'uy', 'Ex', 'Ey'], 'formats': [numpy.float32]*6}
        trajinCI.implicit_dict = {'names': ['x', 'ux', 'phi'], 'formats': [numpy.float32]*3}
        trajinCI.neutral_dict = {'names': ['x', 'ux', 'y', 'uy'], 'formats': [numpy.float32]*4}

        # Add a trajCI reference to the particle object
        pCI = particleCI # abbreviation
        pCI.trajCI = Trajectory_C(trajinCI, ctrlCI, pCI.explicit_species, pCI.implicit_species, pCI.neutral_species)

        ##  Mesh input for the particle mesh, including particle boundary conditions.

        # Create a 2D Cartesian mesh to use for advancing the particles.  The particles
        # themselves are given 3D coordinates.

        from UserMesh_FE_XYZ_Module import UserMeshInput_C

        # 2D mesh input

        umi2DCI = UserMeshInput_C()
        (xmin, ymin) = (-10.0, -10.0)
        (xmax, ymax) = ( 10.0,  10.0)
        umi2DCI.pmin = df_M.Point(xmin, ymin)
        umi2DCI.pmax = df_M.Point(xmax, ymax)
        umi2DCI.cells_on_side = (4, 2)
#        umi2DCI.diagonal = 'crossed'

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

        umi2DCI.particle_boundary_dict = particleBoundaryDict

        ## Create the 2D Cartesian mesh

        from UserMesh_FE_XYZ_Module import UserMesh_C

        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name + ": XY mesh"
        pmeshCI = UserMesh_C(umi2DCI, compute_dictionaries=True, compute_tree=True, plot_flag=False, plot_title=plotTitle)
        # Add this to the particle object:
        particleCI.pmeshCI = pmeshCI


        ### Input for initial particles (i.e., particles present at t=0)

        # a. Name the species (it should be in species_names above)
        speciesName = 'neutral_H'

        # Check that this species has been defined above
        if speciesName not in particleCI.species_names:
            print "The species", speciesName, "has not been defined"
            sys.exit()

        initialDistributionType = 'listed'
        # Check that there's a function listing the particles particles
        printFlag = True
        if hasattr(userParticleClass, speciesName):
            if printFlag: print fncName + "(DnT INFO) Initial distribution for", speciesName, "is the function of that name in", userParticleClass
        # Write error message and exit if no distribution function exists
        else:
            errorMsg = fncName + "(DnT ERROR) Need to define a particle distribution function %s in UserParticle.py for species %s " % (speciesName, speciesName)
            sys.exit(errorMsg)

        # Collect the parameters into a dictionary
        # The 'listed' type will expect a function with the same name as the species.
        neutralHParams = {'species_name': speciesName,
                          'initial_distribution_type': initialDistributionType,
                          }

        # The dictionary keys are mnemonics for initialized particle distributions
        initialParticlesDict = {'initial_neutral_H': (neutralHParams,),
                               }


        # Add the initialized particles to the Particle_C object
        particleCI.initial_particles_dict = initialParticlesDict

# Particle boundary-conditions

        # UserParticleBoundaryFunctions_C is where the facet-crossing callback
        # functions are defined.
#        user_particle_bcs_class = UPrt_M.UserParticleMeshBoundaryConditions_C
        userPBndFnsClass = UPrt_M.UserParticleBoundaryFunctions_C # abbreviation

        # Make the particle-mesh boundary-conditions object and add it
        # to the particle object.  The user has to supply the
        # facet-crossing callback functions in the
        # UserParticleBoundaryFunctions_C object above.

        spNames = particleCI.species_names
        pmeshBCCI = ParticleMeshBoundaryConditions_C(spNames, pmeshCI, userPBndFnsClass, print_flag=False)
        particleCI.pmesh_bcCI = pmeshBCCI

        # Create the initial particles
        printFlags = {}
        for sp in particleCI.species_names: printFlags[sp] = True
        particleCI.initialize_particles(printFlags)

# Get the initial cell index of each particle.

# Should this be something the pmesh computes?  No: pmesh computes the
# index of a single particle.  It doesn't know the particle storage
# infrastructure.
        particleCI.compute_mesh_cell_indices()

# Advance the particles for n_timesteps

        ctrlCI.time_step = 0
        ctrlCI.time = 0.0

        print "Moving", pCI.get_total_particle_count(), "particles for", ctrlCI.n_timesteps, "timesteps"
        for istep in xrange(ctrlCI.n_timesteps):
            particleCI.move_neutral_particles(ctrlCI.dt)

            ctrlCI.time_step += 1
            ctrlCI.time += ctrlCI.dt

        return
#    def test_1_2D_x_y_absorbing_boundary(self):ENDDEF

#class TestParticleBoundaryConditions(unittest.TestCase):
    def test_2_2D_r_theta_absorbing_boundary(self):
        """ Check that particles are deleted correctly when they
            strike an absorbing boundary.
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest: ', fncName

        # if os.environ.get('DISPLAY') is None:
        #     plotFlag=False
        # else:
        #     plotFlag=True

        ### Particle species input

        # Create an instance of the DTparticleInput class
        pinCI = ParticleInput_C()
        # Initialize particles
        pinCI.precision = numpy.float64
        pinCI.particle_integration_loop = 'loop-on-particles'
        pinCI.position_coordinates = ['x', 'y',] # determines the particle-storage dimensions
        pinCI.force_components = ['x', 'y',]
        pinCI.force_precision = numpy.float64

        # Specify the particle-species properties
        # 1. electrons
#         pinCI.particle_species = (('trajelectrons',
#                              {'initial_distribution_type' : 'listed',
#                               'charge' : -1.0*MyPlasmaUnits_C.elem_charge,
#                               'mass' : 1.0*MyPlasmaUnits_C.electron_mass,
#                               'dynamics' : 'explicit',
# #                              'dynamics' : 'implicit',
#                               }
#                              ),
#                             )

        # Define an electron species. Use this name later to initialize this species
        # or to create a source of this species.
        speciesName = 'trajelectrons'
        charge = -1.0*MyPlasmaUnits_C.elem_charge
        mass = 1.0*MyPlasmaUnits_C.electron_mass
        dynamics = 'explicit'
        trajelectronsCI = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        # Add the electrons to particle input
        pinCI.particle_species = (trajelectronsCI,
                                 )
        ## Make the particle storage array for all species.
        particleCI = Particle_C(pinCI, print_flag=False)

        ## Give the name of the .py file containing special particle data (lists of
        # particles, boundary conditions, source regions, etc.)
        userParticleModule = "UserParticles_2D_e"

        # Import this module
        UPrt_M = im_M.import_module(userParticleModule)

        particleCI.user_particle_module = userParticleModule
        particleCI.user_particle_class = userParticleClass = UPrt_M.UserParticleDistributions_C


        ### Add a ref to a Trajectory_C object to particleCI

        # Create input object for trajectories
        trajinCI = TrajectoryInput_C()

        trajinCI.maxpoints = None # Set to None to get every point

        # Specify which particle variables to save.  This has the
        # form of a numpy dtype specification.
        trajinCI.explicit_dict = {'names': ['x', 'ux', 'y', 'uy', 'Ex', 'Ey'], 'formats': [numpy.float32]*6}
        trajinCI.implicit_dict = {'names': ['x', 'ux', 'phi'], 'formats': [numpy.float32]*3}
        trajinCI.neutral_dict = {'names': ['x', 'ux', 'y', 'uy'], 'formats': [numpy.float32]*4}

        # # Initialize the particles
        # printFlags = {}
        # for sp in particleCI.species_names: printFlags[sp] = False
        # particleCI.initialize_particles(printFlags)

        ###  Mesh and Fields input for the particle mesh.

        ## The mesh input
        from UserMesh_y_Fields_FE2D_Module import UserMeshInput_C

        ## The mesh to be created is in an existing file

        umiCI = UserMeshInput_C()
        umiCI.mesh_file = 'quarter_circle_mesh_crossed.xml'
        umiCI.particle_boundary_file='Pbcs_quarter_circle_mesh_crossed.xml'
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

        umiCI.particle_boundary_dict = particleBoundaryDict

        ## Create the particle mesh
        from UserMesh_y_Fields_FE2D_Module import UserMesh_C
        pmeshCI = UserMesh_C(umiCI, compute_dictionaries=True, compute_tree=True, plot_flag=False)

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


        ### Input for initial particles (i.e., particles present at t=0)

        # Name the species (it should be in species_names above)
        speciesName = 'trajelectrons'

        # Check that this species has been defined above
        if speciesName not in particleCI.species_names:
            print "The species", speciesName, "has not been defined"
            sys.exit()

        initialDistributionType = 'listed'
        # Check that there's a function listing the particles particles
        printFlag = True
        if hasattr(userParticleClass, speciesName):
            if printFlag: print fncName + "(DnT INFO) Initial distribution for", speciesName, "is the function of that name in", userParticleClass
        # Write error message and exit if no distribution function exists
        else:
            errorMsg = fncName + "(DnT ERROR) Need to define a particle distribution function %s in UserParticle.py for species %s " % (speciesName, speciesName)
            sys.exit(errorMsg)

        # Collect the parameters into a dictionary
        # The 'listed' type will expect a function with the same name as the species.
        trajElectronParams = {'species_name': speciesName,
                              'initial_distribution_type': initialDistributionType,
                              }

        # The dictionary keys are mnemonics for the initialized particles
        initialParticlesDict = {'initial_trajelectrons': (trajElectronParams,),
                               }

        # Add the initialized particles to the Particle_C object
        particleCI.initial_particles_dict = initialParticlesDict


        ## Particle boundary-conditions

        # UserParticleBoundaryFunctions_C is where the facet-crossing callback
        # functions are defined.
        userPBndFnsClass = UPrt_M.UserParticleBoundaryFunctions_C # abbreviation

        # Make the particle-mesh boundary-conditions object and add it
        # to the particle object.
        spNames = particleCI.species_names
        pmeshBCCI = ParticleMeshBoundaryConditions_C(spNames, pmeshCI, userPBndFnsClass, print_flag=False)
        particleCI.pmesh_bcCI = pmeshBCCI

        ## Set control variables

        ctrlCI = DTcontrol_C()

        # These are fast electrons, so the timestep is small
        ctrlCI.dt = 1.0e-6
        ctrlCI.n_timesteps = 14

        # The trajectory object can now be created and added to particleCI
        pCI = particleCI
        pCI.trajCI = Trajectory_C(trajinCI, ctrlCI, pCI.explicit_species, pCI.implicit_species, pCI.neutral_species)


        # Initialize the particles
        printFlags = {}
        for sp in pCI.species_names: printFlags[sp] = False
        pCI.initialize_particles(printFlags)

        # Get the initial cell index of each particle.
        pCI.compute_mesh_cell_indices()

        ### Particle loop

        print "Moving", pCI.get_total_particle_count(), "particles for", ctrlCI.n_timesteps, "timesteps"

        ctrlCI.time_step = 0
        ctrlCI.time = 0.0
        for istep in xrange(ctrlCI.n_timesteps):

            if pCI.trajCI is not None:
#                print 'pCI.trajCI.skip:', pCI.trajCI.skip
                if istep % pCI.trajCI.skip == 0:
                    pCI.record_trajectory_data(neg_E_field=negElectricField)

            pCI.move_particles_in_electrostatic_field(ctrlCI.dt, negElectricField)

            ctrlCI.time_step += 1
            ctrlCI.time += ctrlCI.dt

        # Record the LAST point on the particle trajectory
        if pCI.trajCI is not None:
                pCI.record_trajectory_data(neg_E_field=negElectricField)

        # Plot the trajectory onto the particle mesh

        mesh = pCI.pmeshCI.mesh
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name + ": XY mesh"
        holdPlot = True # Set to True to stop the plot from disappearing.
        pCI.trajCI.plot_trajectories_on_mesh(mesh, plotTitle, hold_plot=holdPlot) # Plots trajectory spatial coordinates on top of the particle mesh

        return
#    def test_2D_r_theta_absorbing_boundary(self):

#class TestParticleBoundaryConditions(unittest.TestCase):ENDCLASS

if __name__ == '__main__':
    unittest.main()
