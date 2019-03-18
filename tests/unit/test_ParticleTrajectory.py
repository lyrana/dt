#!/usr/bin/env python3

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import os
import numpy
import importlib as im_m
import unittest

import dolfin as df_m

from DT_Module import DTcontrol_C

#from Dolfin_Module import Mesh_C
from Dolfin_Module import Field_C

from Particle_Module import *

from RecordedData_Module import *

from SegmentedArrayPair_Module import SegmentedArrayPair_C

from UserMesh_y_Fields_FE2D_Module import *

from UserUnits_Module import MyPlasmaUnits_C

#STARTCLASS
class TestParticleTrajectory(unittest.TestCase):
    """Test classes in RecordedData_Module.
    """
    
    def setUp(self):

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

        # initializations for each test go here...

        self.ctrl = DTcontrol_C()

        self.ctrl.dt = 1.0e-6

        # Initialize time counters
        self.ctrl.timeloop_count = 0
        self.ctrl.time = 0.0

        # A 'None' value means that the field IS APPLIED
        self.ctrl.apply_solved_electric_field = None
        
        ### Particle species input

        # Create an instance of the DTparticleInput class
        pin = ParticleInput_C()
        # Initialize particles
        pin.precision = numpy.float64
        pin.particle_integration_loop = 'loop-on-particles'
        pin.position_coordinates = ['x', 'y',] # determines the particle-storage dimensions
        pin.force_components = ['x', 'y',]
        pin.force_precision = numpy.float64

        # Define an electron species called 'trajelectrons'.
        speciesName = 'trajelectrons'
        charge = -1.0*MyPlasmaUnits_C.elem_charge
        mass = 1.0*MyPlasmaUnits_C.electron_mass
        dynamics = 'explicit'
        trajelectrons_S = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        # Define an electron species called 'test_electrons'.
        speciesName = 'test_electrons'
        charge = -1.0*MyPlasmaUnits_C.elem_charge
        mass = 1.0*MyPlasmaUnits_C.electron_mass
        dynamics = 'explicit'
        test_electrons_S = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        # Add the electrons to particle input
        pin.particle_species = (trajelectrons_S,
                                test_electrons_S,)

        ## Make the particle storage array for all species.
        self.particle_P = Particle_C(pin, print_flag=True)

        ## Give the name of the .py file containing special particle data (lists of
        # particles, boundary conditions, source regions, etc.)
        userParticlesModuleName = "UserParticles_2D_e"

        # Import this module
        userParticlesModule = im_m.import_module(userParticlesModuleName)

        # self.particle_P.user_particles_module_name = userParticlesModuleName
        self.particle_P.user_particles_class = userParticlesClass = userParticlesModule.UserParticleDistributions_C


        ###  Mesh creation

        umi = UserMeshInput2DCirc_C()

        # Make the mesh & fields from saved files
        umi.mesh_file = 'mesh_quarter_circle_crossed.xml'
        umi.particle_boundary_file='mesh_quarter_circle_crossed_Pbcs.xml'

        ### Input for initial particles (i.e., particles present at t=0)

        # A. trajelectrons
        # Name the species (it should be in species_names above)
        speciesName = 'trajelectrons'

        # Check that this species has been defined above
        if speciesName in self.particle_P.species_names:
            charge = self.particle_P.charge[speciesName]
            mass = self.particle_P.mass[speciesName]
        else:
            print("The species", speciesName, "has not been defined")
            sys.exit()

        initialDistributionType = 'listed'
        # Check that there's a function listing the particles particles
        printFlag = True
        if hasattr(userParticlesClass, speciesName):
            if printFlag: print(fncName + "DnT INFO: Initial distribution for", speciesName, "is the function of that name in", userParticlesClass)
        # Write error message and exit if no distribution function exists
        else:
            errorMsg = fncName + "(DnT ERROR) Need to define a particle distribution function %s in UserParticle.py for species %s " % (speciesName, speciesName)
            sys.exit(errorMsg)

        # Collect the parameters into a dictionary
        # The 'listed' type will expect a function with the same name as the species.
        trajElectronParams = {'species_name': speciesName,
                              'initial_distribution_type': initialDistributionType,
                              }

        # B. test_electrons
        # Name the species (it should be in species_names above)
        speciesName = 'test_electrons'

        # Check that this species has been defined above
        if speciesName in self.particle_P.species_names:
            charge = self.particle_P.charge[speciesName]
            mass = self.particle_P.mass[speciesName]
        else:
            print("The species", speciesName, "has not been defined")
            sys.exit()

        initialDistributionType = 'listed'
        # Check that there's a function listing the particles particles
        printFlag = True
        if hasattr(userParticlesClass, speciesName):
            if printFlag: print(fncName + "DnT INFO: Initial distribution for", speciesName, "is the function of that name in", userParticlesClass)
        # Write error message and exit if no distribution function exists
        else:
            errorMsg = fncName + "(DnT ERROR) Need to define a particle distribution function %s in UserParticle.py for species %s " % (speciesName, speciesName)
            sys.exit(errorMsg)

        # Collect the parameters into a dictionary
        # The 'listed' type will expect a function with the same name as the species.
        test_electronParams = {'species_name': speciesName,
                              'initial_distribution_type': initialDistributionType,
                              }

        ## Collect the initial particles into a dictionary
        
        # The dictionary keys are mnemonics for the initialized particles
        initialParticlesDict = {'initial_trajelectrons': (trajElectronParams,),
                                'initial_test_electrons': (test_electronParams,),
                               }

        # Add the initialized particles to the Particle_C object
        self.particle_P.initial_particles_dict = initialParticlesDict

        # Create the particle mesh object
# 'crossed' diagonals
#        self.pmesh2D = Mesh_C(meshFile="quarter_circle_mesh_crossed.xml", computeDictionaries=True, computeTree=True, plotFlag=False)
# 'left' diagonal
#        self.pmesh2D = Mesh_C(meshFile="quarter_circle_mesh_left.xml", computeDictionaries=True, computeTree=True, plotFlag=False)

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

        # Read the mesh from an existing file

#        pin.pmesh_M = UserMesh_C(meshFile='quarter_circle_mesh_crossed.xml', particleBoundaryFile='Pbcs_quarter_circle_mesh_crossed.xml', computeDictionaries=True, computeTree=True, plotFlag=False)

# Can this be attached to Particle_C after Particle_C construction? YES
        pmesh_M = UserMesh2DCirc_C(umi, compute_dictionaries=True, compute_tree=True, plot_flag=False)

        self.particle_P.pmesh_M = pmesh_M

        ### Field creation

        # The following value should correspond to the element degree
        # used in the potential from which negE was obtained
        phi_element_degree = 1

        if phi_element_degree == 1:
            # For linear elements, grad(phi) is discontinuous across
            # elements. To represent this field, we need Discontinuous Galerkin
            # elements.
            electric_field_element_type = "DG"
        else:
            electric_field_element_type = "Lagrange"

        ## Create the electric field from a file written by test_FieldSolve.py.

        # Put the negative electric field directly on the particle mesh
        self.neg_electric_field = Field_C(pmesh_M,
                                          element_type=electric_field_element_type,
                                          element_degree=phi_element_degree-1,
                                          field_type='vector')

        file = df_m.File("negE2D_crossed-2016.xml") # Use a frozen version.
        file >> self.neg_electric_field.function

        ## Particle boundary-conditions

        # UserParticleBoundaryFunctions_C is where the facet-crossing callback
        # functions are defined.
#        userPBndFnsClass = userParticlesModule.UserParticleBoundaryFunctions_C # LHS is an abbreviation
        userPBndFns = userParticlesModule.UserParticleBoundaryFunctions_C(self.particle_P.position_coordinates, self.particle_P.dx)


        # Make the particle-mesh boundary-conditions object and add it
        # to the particle object.
        spNames = self.particle_P.species_names
#        pmeshBCs = ParticleMeshBoundaryConditions_C(spNames, pmesh_M, userPBndFnsClass, print_flag=False)
        pmeshBCs = ParticleMeshBoundaryConditions_C(spNames, pmesh_M, userPBndFns, print_flag=False)
        self.particle_P.pmesh_bcs = pmeshBCs

        ### Create input for a particle trajectory object

        # Use an input object to collect initialization data for the trajectory object
        self.trajin = TrajectoryInput_C()

        self.trajin.maxpoints = None # Set to None to get every point
        self.trajin.extra_points = 10 # Set to 1 to make sure one boundary-crossing can be
                                      # accommodated. Set to a larger value if there are
                                      # multiple boundary reflections. 10 is needed for
                                      # all the the reflections in test_4.

        # Specify which particle variables to save.  This has the form of a numpy
        # dtype specification.

        name_list_base = ['step', 't'] # These are always recorded
        format_list_base = [int, numpy.float32] # Start off the format list with types for
                                                # 'step' and 't'

        explicit_attributes = ['x', 'ux', 'y', 'uy', 'crossings', 'Ex', 'Ey']
        format_list = format_list_base + [numpy.float32 for i in range(len(explicit_attributes))]
        self.trajin.explicit_dict = {'names': name_list_base+explicit_attributes, 'formats': format_list}
        
        implicit_attributes = ['x', 'ux', 'phi']
        format_list = format_list_base + [numpy.float32 for i in range(len(implicit_attributes))]
        self.trajin.implicit_dict = {'names': name_list_base+implicit_attributes, 'formats': format_list}

        neutral_attributes = ['x', 'ux', 'y', 'uy']
        format_list = format_list_base + [numpy.float32 for i in range(len(neutral_attributes))]
        self.trajin.neutral_dict = {'names': name_list_base+neutral_attributes, 'formats': format_list}

        return

#class TestParticleTrajectory(unittest.TestCase):
    def test_1_trajectory_init(self):
        """ Check that the trajectory variable names are saved correctly in
            the trajectory data_list.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName)

        # Abbreviations
        p_P = self.particle_P

        self.ctrl.n_timesteps = 13

        ## Create the trajectory object and attach it to the particle object.
        # No trajectory storage is created until particles
        # with TRAJECTORY_FLAG on are encountered.
        traj_T = Trajectory_C(self.trajin, self.ctrl, p_P.explicit_species, p_P.implicit_species, p_P.neutral_species)
        self.particle_P.traj_T = traj_T


        # if os.environ.get('DISPLAY') is None:
        #     plotFlag=False
        # else:
        #     plotFlag=True

        # Create particles that are selected for trajectories.
# Could call particle_P.initialize_distributions() instead?

        p_P.create_from_list('trajelectrons', print_flag=True)

        # Check that the trajectory parameters are correct by
        # reading them back from the Particle_C object.

        print("The explicit species with trajectory particles are", p_P.traj_T.explicit_species)
        expectedList = self.trajin.explicit_dict['names']
        for sp in p_P.traj_T.explicit_species:
            print("  with particle indices", p_P.traj_T.particle_index_list[sp])
            if len(p_P.traj_T.particle_index_list[sp]) == 0: continue
            print("  and values", p_P.traj_T.data_list[sp][0].dtype.names)
            gotList = p_P.traj_T.data_list[sp][0].dtype.names
            for i in range(len(expectedList)):
                self.assertEqual(gotList[i], expectedList[i], msg = "Trajectory variable is not correct for an explicit particle")

#        print "The implicit species are", p_P.traj_T.implicit_species
        expectedList = self.trajin.implicit_dict['names']
        for sp in p_P.traj_T.implicit_species:
#            print "  with indices", p_P.traj_T.particle_index_list[sp]
#            print "  and values", p_P.traj_T.data_list[sp][0].dtype.names
            gotList = p_P.traj_T.data_list[sp][0].dtype.names
            for i in range(len(expectedList)):
                self.assertEqual(gotList[i], expectedList[i], msg = "Trajectory variable is not correct for an implicit particle")
#                self.assertAlmostEqual(getparticle[ix], p_expected[isp][ix], msg="Particle is not in correct position")
        return
#    def test_1_trajectory_init(self):ENDDEF

#class TestParticleTrajectory(unittest.TestCase):
    def test_2_record_trajectory(self):
        """ Record and plot the requested trajectory data.
            Check the final particle position.

            The mesh is a 2D quarter circle.
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName)

        p_P = self.particle_P

        self.ctrl.n_timesteps = 13
        dt = self.ctrl.dt

        ## Create the trajectory object and attach it to the particle object.
        # No trajectory storage is created until particles
        # with TRAJECTORY_FLAG on are encountered.
        traj_T = Trajectory_C(self.trajin, self.ctrl, p_P.explicit_species, p_P.implicit_species, p_P.neutral_species)
        self.particle_P.traj_T = traj_T

        # Create particles that are selected for trajectories.
        p_P.create_from_list('trajelectrons', print_flag=True)

        # Get the initial cell index of each particle.
        p_P.compute_mesh_cell_indices()

        # Record the first point on trajectories of marked particles
        if p_P.traj_T is not None:
            p_P.record_trajectory_data(self.ctrl.timeloop_count, self.ctrl.time, neg_E_field=self.neg_electric_field)

        print("\nDnT INFO: %s\t Run for %d timesteps from step %d, time %.3g with %d particles\n" % (fncName, self.ctrl.n_timesteps, self.ctrl.timeloop_count, self.ctrl.time, p_P.get_total_particle_count()))

        for istep in range(self.ctrl.n_timesteps):

            # Count the number of times this loop has been entered:
            self.ctrl.timeloop_count += 1
            # Set the time we're integrating to:
            self.ctrl.time += self.ctrl.dt

            print(fncName, "Starting iteration", self.ctrl.timeloop_count, "to reach time", self.ctrl.time)

            # Do the implicit species first
            if len(p_P.implicit_species) != 0:
                self.iterate_implicit_electrostatic_particles(dt, p_P.implicit_species)

            # Then move the explicit species
            if len(p_P.explicit_species) != 0:
                p_P.move_particles_in_electrostatic_field(self.ctrl, neg_E_field=self.neg_electric_field)

            # XX needs dt; doesn't need n_timesteps

            # Record particle trajectory data
            if p_P.traj_T is not None:
#                print 'p_P.traj_T.skip:', p_P.traj_T.skip
                if self.ctrl.timeloop_count % p_P.traj_T.skip == 0:
                    p_P.record_trajectory_data(self.ctrl.timeloop_count, self.ctrl.time, neg_E_field=self.neg_electric_field)

# test for neg_E_field = None:
#                    p_P.record_trajectory_data()

        print("\nDnT INFO: %s\t Exited the time loop at end of step %d, time %.3g with %d particles\n" % (fncName, self.ctrl.timeloop_count, self.ctrl.time, p_P.get_total_particle_count()))

        # Record the LAST points on the particle trajectories
        if p_P.traj_T is not None:
            p_P.record_trajectory_data(self.ctrl.timeloop_count, self.ctrl.time, neg_E_field=self.neg_electric_field)

        # Check the results

#        print 'trajelectrons data x:', p_P.traj_T.data_list['trajelectrons'][0][:]['x']
#        print 'trajelectrons data ux:', p_P.traj_T.data_list['trajelectrons'][0][:]['ux']
#        print 'trajelectrons data Ex:', p_P.traj_T.data_list['trajelectrons'][0][:]['Ex']

        # First electron
        xp1 = 0.77759792; yp1 = 0.78651935; zp1 = 0.0
        p1 = (xp1,yp1,)

        # Second electron
        xp2 = 1.02059224; yp2 = 0.42618424; zp2 = 0.0
        p2 = (xp2,yp2,)

        p_expected = (p1, p2)

        # Check the results
        ncoords = p_P.particle_dimension # number of particle coordinates to check
        isp = 0
        for sp in p_P.species_names:
#            print 'p_P.get_species_particle_count(sp)', p_P.get_species_particle_count(sp)
            if p_P.get_species_particle_count(sp) == 0: continue

            # Check that the first two particles in the array reaches the correct values
            for ip in [0, 1]:
                getparticle = p_P.pseg_arr[sp].get(ip)
#                print 'calculated = ', getparticle
#                print 'expected = ', p_expected[ip]
                for ic in range(ncoords):
#                    print "result:", getparticle[ic]/p_expected[ip][ic]
# Note: for different field solver, may have to reduce the places:
                    self.assertAlmostEqual(getparticle[ic]/p_expected[ip][ic], 1.0, places=6, msg="Particle is not in correct position")
#                    print "ic", ic, "is OK"
            isp += 1

        # Plot the trajectory onto the particle mesh

        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name
        mesh = p_P.pmesh_M.mesh
        holdPlot = False # Set to True to stop the plot from disappearing.
        p_P.traj_T.plot_trajectories_on_mesh(mesh, plotTitle, hold_plot=holdPlot) # Plots trajectory spatial coordinates on top of the particle mesh

        # Plot the trajectory in phase-space

        plotPhaseSpace = False
        if os.environ.get('DISPLAY') is not None and plotPhaseSpace is True:
            p_P.traj_T.plot() # Phase-space plot of trajectory

        return
#    def test_2_record_trajectory(self):ENDDEF

#class TestParticleTrajectory(unittest.TestCase):
    def test_3_out_of_bounds(self):
        """ Record and plot the requested trajectory data.
            Make the particles leave the mesh.
        """

        printInfoAdvance = True
        printInfoStarting = True
        printInfoExited = True
        printWarningAtEnd = True
        
        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName)

        p_P = self.particle_P

        self.ctrl.n_timesteps = 20 # The particle goes out-of-bounds on step 14
        dt = self.ctrl.dt

        ## Create the trajectory object and attach it to the particle object.
        # No trajectory storage is created until particles
        # with TRAJECTORY_FLAG on are encountered.
        traj_T = Trajectory_C(self.trajin, self.ctrl, p_P.explicit_species, p_P.implicit_species, p_P.neutral_species)
        self.particle_P.traj_T = traj_T

        # Create particles that are selected for trajectories.
        p_P.create_from_list('trajelectrons', print_flag=True)

        # Get the initial cell index of each particle.
        p_P.compute_mesh_cell_indices()

        # Record the first point on trajectories of marked particles
        if p_P.traj_T is not None:
            p_P.record_trajectory_data(self.ctrl.timeloop_count, self.ctrl.time, neg_E_field=self.neg_electric_field)

        if printInfoAdvance is True:
            print("\nDnT INFO: %s\t Advance for %d timesteps from step %d, time %.3g with %d particles\n" % (fncName, self.ctrl.n_timesteps, self.ctrl.timeloop_count, self.ctrl.time, p_P.get_total_particle_count()))

        for istep in range(self.ctrl.n_timesteps):

            self.ctrl.timeloop_count += 1
            self.ctrl.time += self.ctrl.dt

            if printInfoStarting is True:
                print("\nDnT INFO: %s\t Starting iteration %d to reach time %.3g" % (fncName, self.ctrl.timeloop_count, self.ctrl.time))

            # Do the implicit species first
            if len(p_P.implicit_species) != 0:
                self.iterate_implicit_electrostatic_particles(dt, p_P.implicit_species)

            # Then move the explicit species
            if len(p_P.explicit_species) != 0:
                p_P.move_particles_in_electrostatic_field(self.ctrl, neg_E_field=self.neg_electric_field)

            # Gather particle trajectory data for marked particles
            if p_P.traj_T is not None:
#                print 'p_P.traj_T.skip:', p_P.traj_T.skip
                if self.ctrl.timeloop_count % p_P.traj_T.skip == 0:
                    p_P.record_trajectory_data(self.ctrl.timeloop_count, self.ctrl.time, neg_E_field=self.neg_electric_field)

            if p_P.get_total_particle_count() == 0:
                if printWarningAtEnd is True:
                    print("\nDnT WARNING: %s\t At end of step %d, time %.3g, there are 0 particles" % (fncName, self.ctrl.timeloop_count, self.ctrl.time))

        # The specified number of time steps has been reached.
        if printInfoExited is True:
            print("\nDnT INFO: %s\t Exited the time loop at end of step %d, time %.3g with %d particles\n" % (fncName, self.ctrl.timeloop_count, self.ctrl.time, p_P.get_total_particle_count()))

        # Record the LAST points on the particle trajectory
        if p_P.traj_T is not None:
            p_P.record_trajectory_data(self.ctrl.timeloop_count, self.ctrl.time, neg_E_field=self.neg_electric_field)

        # Check the results

#        print 'trajelectron 0 data x:', p_P.traj_T.data_list['trajelectrons'][0][:]['x']
#        print 'trajelectron 1 data x:', p_P.traj_T.data_list['trajelectrons'][1][:]['x']
#        print 'trajelectron 0 data ux:', p_P.traj_T.data_list['trajelectrons'][0][:]['ux']
#        print 'trajelectron 1 data Ex:', p_P.traj_T.data_list['trajelectrons'][0][:]['Ex']

        # The expected results from ParticleNonuniformE.ods
        		
        # First electron
        xp1 = 0.77759792; yp1 = 0.78651935; zp1 = 0.0
        p1 = (xp1,yp1,)

        # Second electron
        xp2 = 1.02059224; yp2 = 0.42618424; zp2 = 0.0
        p2 = (xp2,yp2,)

        p_expected = (p1, p2)

        # Check the results
        ncoords = p_P.particle_dimension # number of particle coordinates to check
        isp = 0
        for sp in p_P.species_names:
            if p_P.get_species_particle_count(sp) == 0:
                print("DnT INFO: %s\t There are no particles remaining for species %s. Skipping check of results." % (fncName, sp))
                continue
            # Check that the first two particles in the array reaches the correct values
            for ip in [0, 1]:
                getparticle = p_P.pseg_arr[sp].get(ip)
#                print 'calculated = ', getparticle
#                print 'expected = ', p_expected[ip]
                for ic in range(ncoords):
#                    print "result:", getparticle[ic]/p_expected[ip][ic]
# Note: for different field solver, may have to reduce the places:
                    self.assertAlmostEqual(getparticle[ic]/p_expected[ip][ic], 1.0, places=6, msg="Particle is not in correct position")
#                    print "ic", ic, "is OK"
            isp += 1

        # Plot the trajectory onto the particle mesh

        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name
        mesh = p_P.pmesh_M.mesh
        holdPlot = True # Set to True to stop the plot from disappearing.
        p_P.traj_T.plot_trajectories_on_mesh(mesh, plotTitle, hold_plot=holdPlot) # Plots trajectory spatial coordinates on top of the particle mesh

        # Plot the trajectory in phase-space

        plotPhaseSpace = True
        if os.environ.get('DISPLAY') is not None and plotPhaseSpace is True:
            p_P.traj_T.plot()

        return
#    def test_3_out_of_bounds(self):ENDDEF

#class TestParticleTrajectory(unittest.TestCase):
    def test_4_reflect_at_boundaries(self):
        """ Record and plot the requested trajectory data.
            Make the particles reflect off the inner boundary.
        """

        printInfoAdvance = True
        printInfoStarting = False
        printWarningAtEnd = True

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName)

        p_P = self.particle_P

        self.ctrl.n_timesteps = 20 # 25; 20; The particle hits the inner boundary on step 14.  Hits the outer boundary on step 21.
        dt = self.ctrl.dt

        ## Create the trajectory object and attach it to the particle object.
        # No trajectory storage is created until particles
        # with TRAJECTORY_FLAG on are encountered.
        traj_T = Trajectory_C(self.trajin, self.ctrl, p_P.explicit_species, p_P.implicit_species, p_P.neutral_species)
        self.particle_P.traj_T = traj_T

        # Create particles that are selected for trajectories.
        p_P.create_from_list('test_electrons', print_flag=True)

        # Get the initial cell index of each particle.
        p_P.compute_mesh_cell_indices()

        # Record the first point on trajectories of marked particles
        if p_P.traj_T is not None:
            p_P.record_trajectory_data(self.ctrl.timeloop_count, self.ctrl.time, neg_E_field=self.neg_electric_field)

        if printInfoAdvance is True:
            print("\nDnT INFO: %s\t Advance for %d timesteps from step %d, time %.3g with %d particles\n" % (fncName, self.ctrl.n_timesteps, self.ctrl.timeloop_count, self.ctrl.time, p_P.get_total_particle_count()))

        for istep in range(self.ctrl.n_timesteps):

            self.ctrl.timeloop_count += 1
            self.ctrl.time += self.ctrl.dt

            if printInfoStarting is True:
                print("\nDnT INFO: %s\t Starting iteration %d to reach time %.3g" % (fncName, self.ctrl.timeloop_count, self.ctrl.time))
#            print fncName, "Starting iteration", self.ctrl.timeloop_count, "to reach time", self.ctrl.time

            # Do the implicit species first
            if len(p_P.implicit_species) != 0:
                self.iterate_implicit_electrostatic_particles(dt, p_P.implicit_species)

            # Then move the explicit species
            if len(p_P.explicit_species) != 0:
                p_P.move_particles_in_electrostatic_field(self.ctrl, neg_E_field=self.neg_electric_field)

            # Gather particle trajectory data for marked particles
            if p_P.traj_T is not None:
#                print 'p_P.traj_T.skip:', p_P.traj_T.skip
                if self.ctrl.timeloop_count % p_P.traj_T.skip == 0:
                    p_P.record_trajectory_data(self.ctrl.timeloop_count, self.ctrl.time, neg_E_field=self.neg_electric_field)

            if p_P.get_total_particle_count() == 0:
                if printWarningAtEnd is True:
                    print("\nDnT WARNING: %s\t At end of step %d, time %.3g, there are 0 particles" % (fncName, self.ctrl.timeloop_count, self.ctrl.time))

        # The specified number of time steps has been reached.
        print("\nDnT INFO: %s\t Exited the time loop at end of step %d, time %.3g with %d particles\n" % (fncName, self.ctrl.timeloop_count, self.ctrl.time, p_P.get_total_particle_count()))

        # Record the LAST points on the particle trajectory
        if p_P.traj_T is not None:
            p_P.record_trajectory_data(self.ctrl.timeloop_count, self.ctrl.time, neg_E_field=self.neg_electric_field)

        # Check the results

        # The expected results from ParticleNonuniformE.ods
        
        # First electron
        xp1 = 0.77759792; yp1 = 0.78651935; zp1 = 0.0
        p1 = (xp1,yp1,)

        # Second electron
        xp2 = 1.02059224; yp2 = 0.42618424; zp2 = 0.0
        p2 = (xp2,yp2,)

        p_expected = (p1, p2)

        # Check the results
        ncoords = p_P.particle_dimension # number of particle coordinates to check
        isp = 0
        for sp in p_P.species_names:
            print(fncName, "Checking results for for species", sp)
            
            if p_P.get_species_particle_count(sp) == 0: continue

            # Check that the first two particles in the array reaches the correct values
            for ip in [0, 1]:
                getparticle = p_P.pseg_arr[sp].get(ip)
#                print 'calculated = ', getparticle
#                print 'expected = ', p_expected[ip]
#                for ic in range(ncoords):
#                    print "result:", getparticle[ic]/p_expected[ip][ic]
# Note: for different field solver, may have to reduce the places:
#                    self.assertAlmostEqual(getparticle[ic]/p_expected[ip][ic], 1.0, places=6, msg="Particle is not in correct position")
#                    print "ic", ic, "is OK"
            isp += 1

        # Plot the trajectory onto the particle mesh

        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name
        mesh = p_P.pmesh_M.mesh
        holdPlot = True # Set to True to stop the plot from disappearing.
        p_P.traj_T.plot_trajectories_on_mesh(mesh, plotTitle, hold_plot=holdPlot) # Plots trajectory spatial coordinates on top of the particle mesh

        # Plot the trajectory in phase-space

        plotPhaseSpace = False
        if os.environ.get('DISPLAY') is not None and plotPhaseSpace is True:
            p_P.traj_T.plot()

        return
#    def test_4_reflect_at_boundaries(self):ENDDEF

#class TestParticleTrajectory(unittest.TestCase):ENDCLASS

if __name__ == '__main__':
    unittest.main()
