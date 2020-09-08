#!/usr/bin/env python3

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
__all__ = ['TestParticleBoundaryConditions.test_1_2D_x_y_absorbing_boundary',
           'TestParticleBoundaryConditions.test_1_cpp_2D_x_y_absorbing_boundary',
           'TestParticleBoundaryConditions.test_2_2D_r_theta_absorbing_boundary',
          ]

import sys
import os
import numpy
import importlib as im_m
import unittest

import dolfin as df_m

from DT_Module import DTcontrol_C

from Dolfin_Module import Field_C

from SegmentedArrayPair_Module import SegmentedArrayPair_C

from Particle_Module import *
from Particle_Module import *

from RecordedData_Module import *

from UserUnits_Module import MyPlasmaUnits_C

#STARTCLASS
class TestParticleBoundaryConditions(unittest.TestCase):
    """Test user-specified boundary conditions on particles.
    """
    
    def setUp(self):

        # Initialization code common to the tests go here...

        self.plot_mesh = False
        self.plot_results = False

        # Turn plots off if there's no display.
        if os.environ.get('DISPLAY') is None:
            self.plot_mesh = False
            self.plot_results = False

        return

#class TestParticleBoundaryConditions(unittest.TestCase):
    def test_1_2D_x_y_absorbing_boundary(self):
        """ Check that particles are deleted correctly when they
            cross an absorbing boundary on a 2D Cartesian mesh.

            The particles are neutral H.
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName)

        if os.environ.get('DISPLAY') is None:
            plotFlag=False
        else:
            plotFlag=self.plot_mesh

        ctrl = DTcontrol_C()

        # Run identifier
        ctrl.title = "test_ParticleBoundaryConditions.py:test_1_2D_x_y_absorbing_boundary"
        # Run author
        ctrl.author = "tph"

        ctrl.timeloop_count = 0
        ctrl.time = 0.0

        ctrl.dt = 0.5
        ctrl.n_timesteps = 100
        ctrl.MAX_FACET_CROSS_COUNT = 100

        ctrl.write_trajectory_files = False

        # Create an instance of the DTparticleInput class
        pin = ParticleInput_C()
        # Settings common to all species
        pin.precision = numpy.float64
        pin.particle_integration_loop = 'loop-on-particles'
        pin.coordinate_system = 'cartesian_xyz'
        pin.force_precision = numpy.float64
        pin.use_cpp_integrators = False # Use C++ code to advance particles.

        # Specify the particle species properties

        speciesName = 'neutral_H' # Use this name later to initialize this
                                  # species or to create a source of this
                                  # species.
        charge = 0.0
        mass = 1.0*MyPlasmaUnits_C.AMU
        dynamics = 'neutral'
        neutralH_S = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        # Add the neutral hydrogen to particle input
        pin.particle_species = (neutralH_S,
                               )

        ## Make the particle storage array for all species.
        particle_P = Particle_C(pin, print_flag=False)

        # Provide the particle distributions for the above species
        # This could be done more like how the mesh is specified:
        # import UserParticles_3D as userParticlesModule
        # Particles have a 3D/3D phase-space even though mesh is just 2D
        userParticlesModuleName = 'UserParticles_3D'

        # Import this module
        infoMsg = "%s\tImporting %s" % (fncName, userParticlesModuleName)
        print(infoMsg)
        userParticlesModule = im_m.import_module(userParticlesModuleName)

        # particle_P.user_particles_module_name = userParticlesModuleName
        particle_P.user_particles_class = userParticlesClass = userParticlesModule.UserParticleDistributions_C

        ### Create a trajectory object and add it to particle_P

        # Create input object for trajectories
        trajin = TrajectoryInput_C()

        trajin.maxpoints = None # Set to None to get every point
        trajin.extra_points = 1 # Set to 1 to make sure one boundary-crossing can be
                                # accommodated.

        # Specify which particle variables to save.  This has the
        # form of a numpy dtype specification.
        name_list_base = ['step', 't'] # These are always recorded
        format_list_base = [int, numpy.float32] # Start off the format list with types for
                                                # 'step' and 't'

        charged_attributes = ['x', 'ux', 'y', 'uy', 'Ex', 'Ey']
        format_list = format_list_base + [numpy.float32]*len(charged_attributes)
        trajin.charged_dict = {'names': name_list_base+charged_attributes, 'formats': format_list}
        
        #implicit_attributes = ['x', 'ux', 'phi']
        #format_list = format_list_base + [numpy.float32]*len(implicit_attributes)
        # trajin.implicit_dict = {'names': name_list_base+implicit_attributes, 'formats': format_list}

        neutral_attributes = ['x', 'ux', 'y', 'uy']
        format_list = format_list_base + [numpy.float32]*len(neutral_attributes)
        trajin.neutral_dict = {'names': name_list_base+neutral_attributes, 'formats': format_list}

        # Add a traj_T reference to the particle object
        p_P = particle_P # abbreviation
        #p_P.traj_T = Trajectory_C(trajin, ctrl, p_P.explicit_species, p_P.implicit_species, p_P.neutral_species)
        p_P.traj_T = Trajectory_C(trajin, ctrl, p_P.charged_species, p_P.neutral_species, p_P.species_index, p_P.mass, p_P.charge)        


        ##  Mesh input for the particle mesh, including particle boundary conditions.

        # Create a 2D Cartesian mesh to use for advancing the particles.  The particles
        # themselves are given 3D coordinates.

        infoMsg = "%s\tImporting UserMeshInput_C from UserMesh_y_Fields_FE_XYZ_Module" % (fncName)
        print(infoMsg)
        from UserMesh_y_Fields_FE_XYZ_Module import UserMeshInput_C
        
        # 2D mesh input

        umi2D = UserMeshInput_C()
        (xmin, ymin) = (-10.0, -10.0)
        (xmax, ymax) = ( 10.0,  10.0)
        umi2D.pmin = df_m.Point(xmin, ymin)
        umi2D.pmax = df_m.Point(xmax, ymax)
        umi2D.cells_on_side = (4, 2)
#        umi2D.diagonal = 'crossed'

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

        umi2D.particle_boundary_dict = particleBoundaryDict

        ## Create the 2D Cartesian mesh
        infoMsg = "%s\tImporting UserMesh_C from UserMesh_y_Fields_FE_XYZ_Module" % (fncName)
        print(infoMsg)
        from UserMesh_y_Fields_FE_XYZ_Module import UserMesh_C

        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name + ": XY mesh"
        pmesh_M = UserMesh_C(umi2D, compute_dictionaries=True, compute_cpp_arrays=False, compute_tree=True, plot_flag=self.plot_mesh, plot_title=plotTitle)

        # Add this to the particle object:
        # p_P.pmesh_M = pmesh_M
        # 1. Attach the particle mesh to p_P.
        # 2. Attach the Python particle movers.
        # 3. Compute the cell-neighbors and facet-normals for the particle movers.
        p_P.initialize_particle_mesh(pmesh_M)
        
        ### Input for initial particles (i.e., particles present at t=0)

        # a. Name the species (it should be in species_names above)
        speciesName = 'neutral_H'

        # Check that this species has been defined above
        if speciesName not in p_P.species_names:
            print("The species", speciesName, "has not been defined")
            sys.exit()

        initialDistributionType = 'listed'
        # Check that there's a function listing the particles particles
        printFlag = True
        if hasattr(userParticlesClass, speciesName):
            if printFlag: print(fncName + "(DnT INFO) Initial distribution for", speciesName, "is the function of that name in", userParticlesClass)
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
        p_P.initial_particles_dict = initialParticlesDict

# Particle boundary-conditions

        # UserParticleBoundaryFunctions_C is where the facet-crossing callback
        # functions are defined by the user.
        # Make an instance of UserParticleBoundaryFunctions_C
        userPBndFns = userParticlesModule.UserParticleBoundaryFunctions_C # abbreviation

        # Make the particle-mesh boundary-conditions object and add it to the
        # particle object.  The user has to supply the facet-crossing callback
        # functions in the UserParticleBoundaryFunctions_C object above.

        spNames = p_P.species_names
        pmeshBCs = ParticleMeshBoundaryConditions_C(spNames, pmesh_M, userPBndFns, print_flag=False)
        p_P.pmesh_bcs = pmeshBCs

        # Create the initial particles
        printFlags = {}
        for sp in p_P.species_names: printFlags[sp] = True
        p_P.initialize_particles(printFlags)

# Get the initial cell index of each particle.

# Should this be something the pmesh computes?  No: pmesh computes the
# index of a single particle.  It doesn't know the particle storage
# infrastructure.
        p_P.compute_mesh_cell_indices()

        p_P.initialize_particle_integration()

# Advance the particles for n_timesteps

        print("Moving", p_P.get_total_particle_count(), "particles for", ctrl.n_timesteps, "timesteps")

        for istep in range(ctrl.n_timesteps):
            
            if p_P.traj_T is not None:
                if istep % p_P.traj_T.skip == 0:
                    p_P.record_trajectory_data(ctrl.timeloop_count, ctrl.time)
                    
            p_P.advance_neutral_particles(ctrl)

            ctrl.timeloop_count += 1
            ctrl.time += ctrl.dt

        # Record the LAST point on the particle trajectory
        if p_P.traj_T is not None:
                p_P.record_trajectory_data(ctrl.timeloop_count, ctrl.time)

        # Plot the trajectory onto the particle mesh
        if self.plot_results is True:
            mesh = p_P.pmesh_M.mesh
            plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name + ": XY mesh"
            holdPlot = True # Set to True to stop the plot from disappearing.
            p_P.traj_T.plot_trajectories_on_mesh(mesh, plotTitle, hold_plot=holdPlot) # Plots trajectory spatial coordinates on top of the particle mesh
            
        return
#    def test_1_2D_x_y_absorbing_boundary(self):ENDDEF

#class TestParticleBoundaryConditions(unittest.TestCase):
    def test_2_cpp_2D_x_y_absorbing_boundary(self):
        """ Check that particles are deleted correctly when they
            cross an absorbing boundary on a 2D Cartesian mesh.

            The particles are neutral H.
            C++ code is used to advance the particles.
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName)

        if os.environ.get('DISPLAY') is None:
            plotFlag=False
        else:
            plotFlag=self.plot_mesh

        ctrl = DTcontrol_C()

        # Run identifier
        ctrl.title = "test_ParticleBoundaryConditions.py:test_2_cpp_2D_x_y_absorbing_boundary"
        # Run author
        ctrl.author = "tph"

        # Set the start time
        ctrl.timeloop_count = 0
        ctrl.time = 0.0

        # Set the timestep and running time
        ctrl.dt = 0.5
        ctrl.n_timesteps = 100 # 100
        ctrl.MAX_FACET_CROSS_COUNT = 100

        ctrl.write_trajectory_files = False

        # Create an instance of the DTparticleInput class
        pin = ParticleInput_C()
        # Settings common to all species
        pin.precision = numpy.float64
        pin.particle_integration_loop = 'loop-on-particles'

        # Set the particle coordinate system. This determines the particle-storage
        # dimensions
        pin.coordinate_system = 'cartesian_xyz'
        pin.force_precision = numpy.float64
        pin.use_cpp_integrators = True # Use C++ code to advance particles.

        # Specify the particle species properties

        speciesName = 'neutral_H' # Use this name later to initialize this
                                  # species and/or to create a source of this
                                  # species.
        charge = 0.0
        mass = 1.0*MyPlasmaUnits_C.AMU
        dynamics = 'neutral'
        neutralH_S = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        # Add the neutral hydrogen to particle input
        pin.particle_species = (neutralH_S,
                               )

        ## Make the particle storage array for all species.
        particle_P = Particle_C(pin, print_flag=False)

        # Provide particle distributions for the above species
        # This could be done more like how the mesh is specified:
        # import UserParticles_3D as userParticlesModule
        
        # Particles have a 3D/3D phase-space even though mesh is just 2D
        userParticlesModuleName = 'UserParticles_3D'

        # Import this module
        userParticlesModule = im_m.import_module(userParticlesModuleName)

        # particle_P.user_particles_module_name = userParticlesModuleName
        particle_P.user_particles_class = userParticlesClass = userParticlesModule.UserParticleDistributions_C

        ### Create a trajectory object and add it to particle_P
        # This will contain all the trajectory data recorded during the simulation.

        # Create input object for trajectories
        trajin = TrajectoryInput_C()

        trajin.maxpoints = None # Set to None to get every point
        trajin.extra_points = 1 # Set to 1 to make sure one boundary-crossing can be
                                # accommodated.

        # Specify which particle variables to save.  This has the form of a numpy
        # dtype specification.
        name_list_base = ['step', 't'] # These variables are always recorded
        format_list_base = [int, numpy.float32] # Start off the format-spec list with
                                                # types for 'step' and 't'

        charged_attributes = ['x', 'ux', 'y', 'uy', 'Ex', 'Ey']
        format_list = format_list_base + [numpy.float32]*len(charged_attributes)
        trajin.charged_dict = {'names': name_list_base+charged_attributes, 'formats': format_list}
        
        #implicit_attributes = ['x', 'ux', 'phi']
        #format_list = format_list_base + [numpy.float32]*len(implicit_attributes)
        # trajin.implicit_dict = {'names': name_list_base+implicit_attributes, 'formats': format_list}

        neutral_attributes = ['x', 'ux', 'y', 'uy']
        format_list = format_list_base + [numpy.float32]*len(neutral_attributes)
        trajin.neutral_dict = {'names': name_list_base+neutral_attributes, 'formats': format_list}

        # Add a traj_T reference to the particle object
        p_P = particle_P # abbreviation
        #p_P.traj_T = Trajectory_C(trajin, ctrl, p_P.explicit_species, p_P.implicit_species, p_P.neutral_species)
        p_P.traj_T = Trajectory_C(trajin, ctrl, p_P.charged_species, p_P.neutral_species, p_P.species_index, p_P.mass, p_P.charge)

        ##  Mesh input for the particle mesh, including particle boundary conditions.

        # Create a 2D Cartesian mesh to use for advancing the particles.  The particles
        # themselves are given 3D coordinates.

        from UserMesh_y_Fields_FE_XYZ_Module import UserMeshInput_C

        # Make input for a square 2D mesh

        umi2D = UserMeshInput_C()
        (xmin, ymin) = (-10.0, -10.0)
        (xmax, ymax) = ( 10.0,  10.0)
        umi2D.pmin = df_m.Point(xmin, ymin)
        umi2D.pmax = df_m.Point(xmax, ymax)
        umi2D.cells_on_side = (4, 2)
#        umi2D.diagonal = 'crossed'

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

        umi2D.particle_boundary_dict = particleBoundaryDict

        ## Create the 2D Cartesian mesh

        from UserMesh_y_Fields_FE_XYZ_Module import UserMesh_C

        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name + ": XY mesh"
        pmesh_M = UserMesh_C(umi2D, compute_dictionaries=True, compute_cpp_arrays=False, compute_tree=True, plot_flag=self.plot_mesh, plot_title=plotTitle)

        # Add this to the particle object:
        # p_P.pmesh_M = pmesh_M
        # 1. Attach the particle mesh to p_P.
        # 2. Attach the C++ particle movers.
        # 3. Compute the cell-neighbors and facet-normals for the particle movers.
        p_P.initialize_particle_mesh(pmesh_M)
        
        ### Input for initial particles (i.e., particles present at t=0)

        # a. Name the species (it has to be one of the species_names above)
        speciesName = 'neutral_H'
        # Check that this species has been defined above
        if speciesName not in p_P.species_names:
            print("The species", speciesName, "has not been defined")
            sys.exit() # change to raiseError

        initialDistributionType = 'listed'
        # Check that there's a function listing the particles particles
        printFlag = True
        if hasattr(userParticlesClass, speciesName):
            if printFlag: print(fncName + "(DnT INFO) Initial distribution for", speciesName, "is the function of that name in", userParticlesClass)
        # Write error message and exit if no distribution function exists
        else:
            errorMsg = fncName + "(DnT ERROR) Need to define a particle distribution function %s in UserParticle.py for species %s " % (speciesName, speciesName)
            raise RuntimeError(errorMsg) # sys.exit(errorMsg)

        # Collect the parameters into a dictionary:
        # A 'listed' type will expect a function with the same name as the species.
        neutralHParams = {'species_name': speciesName,
                          'initial_distribution_type': initialDistributionType,
                         }

        # Collect the initialized species into a dictionary:
        # The dictionary keys are mnemonics for initialized particle distributions
        initialParticlesDict = {'initial_neutral_H': (neutralHParams,),
                               }

        # Add the initialized particles to the Particle_C object
        p_P.initial_particles_dict = initialParticlesDict

        # Add particle boundary-conditions

        # See the TODO of 4apr20 for redoing this.
        
        # Make the particle-mesh boundary-conditions object and add it to the
        # particle object.  The user has to supply the facet-crossing callback
        # functions in the UserParticleBoundaryFunctions_C object above.

        # Import C++ particle boundary-condition functions
        userParticleBoundaryFunctionsSOlibName = "user_particle_boundary_functions_cartesian_xyz_solib"
        infoMsg = "%s\tImporting %s" % (fncName, userParticleBoundaryFunctionsSOlibName)
        print(infoMsg)
        userParticleBoundaryFunctionsSOlib = im_m.import_module(userParticleBoundaryFunctionsSOlibName)
        # Call the constructor to make a UserParticleBoundaryFunctions object
        userPBndFns = userParticleBoundaryFunctionsSOlib.UserParticleBoundaryFunctions(p_P.position_coordinates)
        # Create the map from mesh facets to particle callback functions:
        spNames = p_P.species_names
        pmeshBCs = p_P.particle_solib.ParticleMeshBoundaryConditions(spNames, pmesh_M, userPBndFns, print_flag=False)        

        # Add pmeshBCs to the Particle_C object
        p_P.pmesh_bcs = pmeshBCs
        
        # Create the initial particles
        printFlags = {}
        for sp in p_P.species_names: printFlags[sp] = True
        p_P.initialize_particles(printFlags)

# Get the initial cell index of each particle.

# Should this be something the pmesh computes?  No: pmesh computes the
# index of a single particle.  It doesn't know the particle storage
# infrastructure.
        p_P.compute_mesh_cell_indices()

        p_P.initialize_particle_integration()

# Advance the particles for n_timesteps

        print("Moving", p_P.get_total_particle_count(), "particles for", ctrl.n_timesteps, "timesteps")

        for istep in range(ctrl.n_timesteps):
            
            if p_P.traj_T is not None:
                if istep % p_P.traj_T.skip == 0:
                    p_P.record_trajectory_data(ctrl.timeloop_count, ctrl.time)
                    
            p_P.advance_neutral_particles(ctrl)

            ctrl.timeloop_count += 1
            ctrl.time += ctrl.dt

        # Record the LAST point on the particle trajectory
        if p_P.traj_T is not None:
                p_P.record_trajectory_data(ctrl.timeloop_count, ctrl.time)

        # Plot the trajectory onto the particle mesh
        if self.plot_results is True:
            mesh = p_P.pmesh_M.mesh
            plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name + ": XY mesh"
            holdPlot = True # Set to True to stop the plot from disappearing.
            p_P.traj_T.plot_trajectories_on_mesh(mesh, plotTitle, hold_plot=holdPlot) # Plots trajectory spatial coordinates on top of the particle mesh
            
        return
#    def test_2_cpp_2D_x_y_absorbing_boundary(self):ENDDEF

#class TestParticleBoundaryConditions(unittest.TestCase):
    def test_3_2D_r_theta_absorbing_boundary(self):
        """ Check that particles are deleted correctly when they
            cross a 2D absorbing boundary.

            The particles are electrons and there's an applied Electric field
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName)

        if os.environ.get('DISPLAY') is None:
            plotFlag=False
        else:
            plotFlag=self.plot_results

        ## Set control variables

        ctrl = DTcontrol_C()

        # Run identifier
        ctrl.title = "test_ParticleBoundaryConditions.py:test_3_2D_r_theta_absorbing_boundary"
        # Run author
        ctrl.author = "tph"

        ctrl.timeloop_count = 0
        ctrl.time = 0.0

        # These are fast electrons, so the timestep is small
        ctrl.dt = 1.0e-6
        ctrl.n_timesteps = 14
        ctrl.MAX_FACET_CROSS_COUNT = 100

        ctrl.write_trajectory_files = False

        ### Particle species input

        # Create an instance of the DTparticleInput class
        pin = ParticleInput_C()
        # Initialize particles
        pin.precision = numpy.float64
        pin.particle_integration_loop = 'loop-on-particles'
        pin.coordinate_system = 'cartesian_xy'        
        pin.force_components = ['x', 'y',]
        pin.force_precision = numpy.float64
        pin.use_cpp_integrators = False # Use C++ code to advance particles.
        
        # Specify the particle-species properties

        # Define an electron species. Use this name later to initialize this species
        # or to create a source of this species.
        speciesName = 'trajelectrons'
        charge = -1.0*MyPlasmaUnits_C.elem_charge
        mass = 1.0*MyPlasmaUnits_C.electron_mass
        dynamics = 'explicit'
        trajelectrons_S = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        # Add the electrons to particle input
        pin.particle_species = (trajelectrons_S,
                               )
        ## Make the particle storage array for all species.
        particle_P = Particle_C(pin, print_flag=False)

        ## Give the name of the .py file containing special particle data (lists of
        # particles, boundary conditions, source regions, etc.)
        userParticlesModuleName = "UserParticles_2D_e"

        # Import this module
        userParticlesModule = im_m.import_module(userParticlesModuleName)

        # particle_P.user_particles_module_name = userParticlesModuleName
        particle_P.user_particles_class = userParticlesClass = userParticlesModule.UserParticleDistributions_C


        ### Add a ref to a Trajectory_C object to particle_P

        # Create input object for trajectories
        trajin = TrajectoryInput_C()

        trajin.maxpoints = None # Set to None to get every point
        trajin.extra_points = 1 # Set to 1 to make sure one boundary-crossing can be
                                # accommodated.

        # Specify which particle variables to save ('step' and 't' are required). This has the
        # form of a numpy dtype specification.

        format_list_base = [int]
        format_list = format_list_base + [numpy.float32]*7
        trajin.charged_dict = {'names': ['step', 't', 'x', 'ux', 'y', 'uy', 'Ex', 'Ey'], 'formats': format_list}
        format_list = format_list_base + [numpy.float32]*4
        trajin.implicit_dict = {'names': ['step', 't', 'x', 'ux', 'phi'], 'formats': format_list}
        format_list = format_list_base + [numpy.float32]*5
        trajin.neutral_dict = {'names': ['step', 't', 'x', 'ux', 'y', 'uy'], 'formats': format_list}

        ###  Mesh and Fields input for the particle mesh.

        ## The mesh input
        from UserMesh_y_Fields_FE2D_Module import UserMeshInput2DCirc_C

        ## The mesh to be created is in an existing file

        umi = UserMeshInput2DCirc_C()
        umi.mesh_file = 'mesh_quarter_circle_crossed.xml'
        umi.particle_boundary_file='mesh_quarter_circle_crossed_Pbcs.xml'
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

        umi.particle_boundary_dict = particleBoundaryDict

        ## Create the particle mesh
        from UserMesh_y_Fields_FE2D_Module import UserMesh2DCirc_C
        pmesh_M = UserMesh2DCirc_C(umi, compute_dictionaries=True, compute_cpp_arrays=False, compute_tree=True, plot_flag=self.plot_mesh)

        # Add this to the particle object:
        particle_P.pmesh_M = pmesh_M

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
        negElectricField = Field_C(particle_P.pmesh_M,
                                   element_type=electricFieldElementType,
                                   element_degree=phiElementDegree-1,
                                   field_type='vector')

        file = df_m.File('negE_test_2_2D.xml')
        file >> negElectricField.function


        ### Input for initial particles (i.e., particles present at t=0)

        # Name the species (it should be in species_names above)
        speciesName = 'trajelectrons'

        # Check that this species has been defined above
        if speciesName not in particle_P.species_names:
            print("The species", speciesName, "has not been defined")
            sys.exit()

        initialDistributionType = 'listed'
        # Check that there's a function listing the particles particles
        printFlag = True
        if hasattr(userParticlesClass, speciesName):
            if printFlag: print(fncName + "(DnT INFO) Initial distribution for", speciesName, "is the function of that name in", userParticlesClass)
        # Write error message and exit if no distribution function exists
        else:
            errorMsg = fncName + "(DnT ERROR) Need to define a particle distribution function %s in UserParticle.py for species %s " % (speciesName, speciesName)
            sys.exit(errorMsg)

        # Collect the parameters into a dictionary

        # For a 'listed' type, there needs to be a function with the same name as the
        # species in userParticlesClass
        trajElectronParams = {'species_name': speciesName,
                              'initial_distribution_type': initialDistributionType,
                              }

        # The dictionary keys are mnemonics strings for identifying each set of
        # initialized particles. They're not the names of particle species.
        initialParticlesDict = {'initial_trajelectrons': (trajElectronParams,),
                               }

        # Add the initialized particles to the Particle_C object
        particle_P.initial_particles_dict = initialParticlesDict


        ## Particle boundary-conditions

        # UserParticleBoundaryFunctions_C is where the facet-crossing callback
        # functions are defined.
        userPBndFnsClass = userParticlesModule.UserParticleBoundaryFunctions_C # abbreviation

        # Make the particle-mesh boundary-conditions object and add it
        # to the particle object.
        spNames = particle_P.species_names
        pmeshBCs = ParticleMeshBoundaryConditions_C(spNames, pmesh_M, userPBndFnsClass, print_flag=False)
        particle_P.pmesh_bcs = pmeshBCs

        # The trajectory object can now be created and added to particle_P
        p_P = particle_P
        # p_P.traj_T = Trajectory_C(trajin, ctrl, p_P.explicit_species, p_P.implicit_species, p_P.neutral_species)
        p_P.traj_T = Trajectory_C(trajin, ctrl, p_P.charged_species, p_P.neutral_species, p_P.species_index, p_P.mass, p_P.charge)        

        # Initialize the particles
        printFlags = {}
        for sp in p_P.species_names: printFlags[sp] = False
        p_P.initialize_particles(printFlags)

        # Get the initial cell index of each particle.
        p_P.compute_mesh_cell_indices()

        p_P.initialize_particle_integration()
        
        ### Particle loop

        print("Moving", p_P.get_total_particle_count(), "particles for", ctrl.n_timesteps, "timesteps")

        for istep in range(ctrl.n_timesteps):

#            print("istep", istep)
            
            if p_P.traj_T is not None:
#                print 'p_P.traj_T.skip:', p_P.traj_T.skip
                if istep % p_P.traj_T.skip == 0:
                    p_P.record_trajectory_data(ctrl.timeloop_count, ctrl.time, neg_E_field=negElectricField)

            p_P.advance_charged_particles_in_E_field(ctrl, neg_E_field=negElectricField)

            ctrl.timeloop_count += 1
            ctrl.time += ctrl.dt

        # Record the LAST point on the particle trajectory
        if p_P.traj_T is not None:
                p_P.record_trajectory_data(ctrl.timeloop_count, ctrl.time, neg_E_field=negElectricField)

        # Plot the trajectory onto the particle mesh
        if self.plot_results is True:
            mesh = p_P.pmesh_M.mesh
            plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name + ": XY mesh"
            holdPlot = True # Set to True to stop the plot from disappearing.
            p_P.traj_T.plot_trajectories_on_mesh(mesh, plotTitle, hold_plot=holdPlot) # Plots trajectory spatial coordinates on top of the particle mesh

        return
#     def test_3_2D_r_theta_absorbing_boundary(self):ENDDEF

#class TestParticleBoundaryConditions(unittest.TestCase):ENDCLASS

if __name__ == '__main__':
    unittest.main()
