#!/usr/bin/env python3

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import os
import numpy as np_m
import importlib as im_m
import unittest

import dolfin as df_m
import matplotlib.pyplot as mplot_m

from DT_Module import DTcontrol_C
from UserUnits_Module import MyPlasmaUnits_C
from Particle_Module import *

# Here's the mesh definition for this test
from UserMesh_y_Fields_FE_XYZ_Module import *

#STARTCLASS
class TestCppParticleMigration(unittest.TestCase):
    """Test use of pybind11 to interface with C++

       Use a C++ function to advance particles, passing needed parameters from Python.
       Neutral helium is the only species.  The initial particles are specified in
       UserParticles_3D.py.

    """
    
    def setUp(self):

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

        # Initializations performed before each test go here...

        self.plot_results = False
        
        # Create an instance of the DTparticleInput class
        pin = ParticleInput_C()

        # Set up particle variables
        pin.precision = np_m.float64
        pin.particle_integration_loop = 'loop-on-particles'

        # The particle coordinate system is 3D Cartesian
        pin.coordinate_system = 'cartesian_xyz'

        # Neutral particles: No electric forces.
        """
        pin.force_components = ['x', 'y',]
        """
        pin.force_precision = np_m.float64
        pin.use_cpp_integrators = True # Use C++ version of particle movers.

        # Give the properties of the particle species.  The charges and masses are
        # normally those of the physical particles, and not the computational
        # macroparticles.  Macroparticle weights are specified or computed in a
        # separate file (see user_particles_module_name below) specifying the
        # particle distribution functions, and can vary from particle to particle.

        speciesName = 'neutral_H'
        charge = 0.0
        mass = 1.0*MyPlasmaUnits_C.AMU
        dynamics = 'neutral'
#        integratorName = "advance_neutral_species"        
#        neutralH_S = ParticleSpecies_C(speciesName, charge, mass, dynamics, integratorName)
        neutralH_S = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        # Add these species to particle input
        pin.particle_species = (neutralH_S,
                               )
        # Make the particle object from pin
        self.particle_P = Particle_C(pin, print_flag=False)

        # Give the name of the .py file containing additional particle data (lists of
        # particles, boundary-condition callbacks, source regions, etc.)
        # Particles have a 3D/3D phase-space even though mesh is just 1D
        userParticlesModuleName = "UserParticles_3D"
#        print "setUp: UserParticle module is", userParticlesModuleName

        # Import this module
        userParticlesModule = im_m.import_module(userParticlesModuleName)

        # self.particle_P.user_particles_module_name = userParticlesModuleName
        self.particle_P.user_particles_class = userParticlesClass = userParticlesModule.UserParticleDistributions_C

        ### Neutral H atoms are present at t=0
        speciesName = 'neutral_H'
        # Check that this species has been defined above
        if speciesName not in self.particle_P.species_names:
            print(fncName + "The species", speciesName, "has not been defined")
            sys.exit()

        # Specify how the species will be initialized
        initialDistributionType = 'listed'
        # Check that there's a function listing the particles particles
        printFlag = False
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

        # The dictionary keys are mnemonics for the initialized particles
        initialParticlesDict = {
                                'initial_neutral_H': (neutralHParams,),
                               }

        self.particle_P.initial_particles_dict = initialParticlesDict

        # Create two initial particles from the above specs.
        for ip in initialParticlesDict:
            ipList = initialParticlesDict[ip]
            ipParams = ipList[0] # The first list item is a dictionary of "particle parameters".
            s = ipParams['species_name'] # Look up the name of the species.
            initialDistributionType = ipParams['initial_distribution_type'] # Look up how the particle distribution will be specified.
            if initialDistributionType == 'listed':
                # Put user-listed particles into the storage array
                self.particle_P.create_from_list(s, False)

        plotFlag = False
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name

        # Turn off plotting if there's no DISPLAY

        # if os.environ.get('DISPLAY') is None:
        #     plotFlag=False
        # else:
        #     plotFlag=True

        # Create a 1D mesh that the particles can move through

        # 1d mesh input
        umi1D = UserMeshInput_C()
        umi1D.pmin = df_m.Point(-10.0)
        umi1D.pmax = df_m.Point(10.0)
        umi1D.cells_on_side = (4,)

        # Identify where particle boundary-conditions will be imposed
        xminIndx = 1
        xmaxIndx = 2
        particleBoundaryDict = {'xmin': xminIndx,
                                'xmax': xmaxIndx,
                                }

        umi1D.particle_boundary_dict = particleBoundaryDict

        # Create a 1D particle mesh
        self.pmesh1D = UserMesh_C(umi1D, compute_dictionaries=False, compute_cpp_arrays=False, compute_tree=True, plot_flag=plotFlag, plot_title=plotTitle + ": 1D")

# Individual geometry dictionaries:
#        self.pmesh1D.compute_cell_vertices_dict()
#        self.pmesh1D.compute_cell_dict()

        # 2D mesh input
        umi2D = UserMeshInput_C()
        umi2D.pmin = df_m.Point(-10.0, -10.0)
        umi2D.pmax = df_m.Point(10.0, 10.0)
        umi2D.cells_on_side = (4, 2)

        # Identify where particle boundary-conditions will be imposed
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

        # Create a 2D particle mesh
        self.pmesh2D = UserMesh_C(umi2D, compute_dictionaries=False, compute_cpp_arrays=False, compute_tree=True, plot_flag=plotFlag, plot_title=plotTitle + ": 2D")
#        self.pmesh2D.compute_cell_vertices_dict()
#        self.pmesh2D.compute_cell_dict()

        # 3D mesh input
        umi3D = UserMeshInput_C()
        umi3D.pmin = df_m.Point(-10.0, -10.0, -10.0)
        umi3D.pmax = df_m.Point(10.0, 10.0, 10.0)
        umi3D.cells_on_side = (4, 4, 4)

        # Identify where particle boundary-conditions will be imposed
        xminIndx = 1
        xmaxIndx = 2
        yminIndx = 4
        ymaxIndx = 8
        zminIndx = 16
        zmaxIndx = 32
        particleBoundaryDict = {'xmin': xminIndx,
                                'xmax': xmaxIndx,
                                'ymin': yminIndx,
                                'ymax': ymaxIndx,
                                'zmin': zminIndx,
                                'zmax': zmaxIndx,
                                }

        umi3D.particle_boundary_dict = particleBoundaryDict


        # Create a 3D particle mesh
        self.pmesh3D = UserMesh_C(umi3D, compute_tree=True, plot_flag=plotFlag, plot_title=plotTitle + ": 3D")
        # Explicitly compute dictionaries needed
#        self.pmesh3D_M.compute_cell_entity_indices_dict('vertex')
#        self.pmesh3D_M.compute_cell_vertices_dict()
#        self.pmesh3D_M.compute_cell_entity_indices_dict('facet')
#        self.pmesh3D_M.compute_cell_dict()
#        self.pmesh3D_M.compute_cell_facet_normals_dict()
#        self.pmesh3D_M.compute_cell_neighbors_dict()

        return
#    def setUp(self):ENDDEF

    def test_1D_particle_migration(self):
        """Test that we can call a C++ function to push neutral particles.

           Initial particle positions are in UserParticles_3D.py.

           The first particle moves in -x only, from 9.5 to -9.5.
           The first particle starts at:
               (x0, y0, z0) = (9.5, -9.5, 0.0), with velocity:
               (ux0, uy0, uz0) = (-2.0, 0.0, 0.0)
           It moves to (-9.5, -9.5, 0.0)

           The second particle starts at
               (x1, y1, z1) = (9.5, 9.5, 9.5)
               (ux1, uy1, uz1) = (-2.0, -2.0, -2.0)
           It moves in -x, -y, -z to the opposite corner.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName)

        ctrl = DTcontrol_C()

        # Run identifier
        ctrl.title = "1D particle advance using C++"
        # Run author
        ctrl.author = "tph"

        ctrl.time = 0.0
        ctrl.timeloop_count = 0

        ctrl.dt = 0.5
        ctrl.n_timesteps = 19
        ctrl.MAX_FACET_CROSS_COUNT = 100

        # Run on a 1D mesh

        # 1. Attach the particle mesh to particle_P
        # 2. Attach the C++ particle movers.
        # 3. Compute the cell-neighbors and facet-normals for the particle movers.
        self.particle_P.initialize_particle_mesh(self.pmesh1D)

        # Attach the particle integrators
        self.particle_P.initialize_particle_integration()

        ### Particle boundary-conditions

        # See the TODO of 4apr20 for redoing this.
        
        # Make a dictionary associating the above-named boundaries of the particle mesh with
        # user-supplied call-back functions.
        
        # First, we need to make a UserParticleBoundaryFunctions_... object, which
        # contains the boundary-crossing callback functions.
        
        spNames = self.particle_P.species_names
        
        # Import C++ particle boundary-conditions
        userParticleBoundaryFunctionsSOlibName = "user_particle_boundary_functions_cartesian_xyz_solib"
        if userParticleBoundaryFunctionsSOlibName not in sys.modules:        
            infoMsg = "%s\t\"Importing %s\"" % (fncName, userParticleBoundaryFunctionsSOlibName)
            print(infoMsg)
        userParticleBoundaryFunctionsSOlib = im_m.import_module(userParticleBoundaryFunctionsSOlibName)
        # Call the constructor to make a UserParticleBoundaryFunctions_... object
        userPBndFns = userParticleBoundaryFunctionsSOlib.UserParticleBoundaryFunctions(self.particle_P.position_coordinates)
        # Create the map from mesh facets to particle callback functions:
        pmeshBCs = self.particle_P.particle_solib.ParticleMeshBoundaryConditions(spNames, self.particle_P.pmesh_M, userPBndFns, print_flag=False)

        # Add pmeshBCs to the Particle_C object
        self.particle_P.pmesh_bcs = pmeshBCs

        # Get the initial cell index of each particle.
        
        # Note: this could be moved down closer to the beginning of the
        # particle-advance loop.
        # cell_dict{} is needed by the Python version of is_inside_cell(), which is
        # used in compute_mesh_cell_indices() below.
        self.particle_P.pmesh_M.compute_cell_dict()
        self.particle_P.compute_mesh_cell_indices()

        ### Put the expected ending results for the two particles into the p_expected tuple ###

        # First particle

        xsp0 = -9.5; ysp0 = -9.5; zsp0 = 0.0
        xsp00 = -8.5; ysp00 = -9.5; zsp00 = 0.0
        vxsp0 = -2.0; vysp0 = 0.0; vzsp0 = 0.0

        weight0 = 2.0
        bitflag0 = 2
        cell_index0 = 0
        unique_ID0 = 0
        crossings = 0

        psp0 = (xsp0,ysp0,zsp0, xsp00,ysp00,zsp00, vxsp0,vysp0,vzsp0, weight0, bitflag0, cell_index0, unique_ID0, crossings)

        # Second particle

        xsp1 = -9.5; ysp1 = -9.5; zsp1 = -9.5
        xsp10 = -8.5; ysp10 = -8.5; zsp10 = -8.5
        vxsp1 = -2.0; vysp1 = -2.0; vzsp1 = -2.0

        weight1 = 3.0
        bitflag1 = 2
        cell_index1 = 0
        unique_ID1 = 1
        crossings = 0

        psp1 = (xsp1,ysp1,zsp1, xsp10,ysp10,zsp10, vxsp1,vysp1,vzsp1, weight1, bitflag1, cell_index1, unique_ID1, crossings)

        p_expected = (psp0, psp1)

        #
        # Move the particles and check the final positions
        #
        
        speciesName = 'neutral_H'
        sap = self.particle_P.sap_dict[speciesName] # segmented array for this species

        # Print the cell indices
#        print("cell_vertices_dict = ", self.particle_P.pmesh_M.cell_vertices_dict)

        # Integrate for n_timesteps
        print("Advancing", self.particle_P.get_total_particle_count(), "particles for", ctrl.n_timesteps, "timesteps on a 1D mesh")
        for istep in range(ctrl.n_timesteps):
#            print(fncName, "istep:", istep)
            self.particle_P.advance_neutral_particles(ctrl)

        # Check the results
        ncoords = self.particle_P.particle_dimension # number of particle coordinates to check
        for sp in self.particle_P.neutral_species:
            for ip in [0, 1]:
                # getparticle = self.particle_P.sap_dict[sp].get(ip)
                # Instead of .get(), retrieve the particle structure using the returned Numpy array it's in.
                (pseg, offset) = self.particle_P.sap_dict[sp].get_segment_and_offset(ip)
                getparticle = pseg[offset] # Retrieve the particle from the SAP.
#                print('expected = ', p_expected[ip])
#                print('calculated = ', getparticle)
                for ic in range(2*ncoords):
                    self.assertAlmostEqual(p_expected[ip][ic], getparticle[ic], places=6, msg="Particle is not in correct position")
                cell_index_position = -3
#                print fncName, "expected cell =", p_expected[ip][cell_index_position], "computed cell =", getparticle[cell_index_position]
                self.assertEqual(p_expected[ip][cell_index_position], getparticle[cell_index_position], msg="Particle is not in correct cell")

        return
    # ENDDEF: def test_1D_particle_migration(self)

    def test_2D_particle_migration(self):
        """Test particle migration across cells on 2D mesh.

        """
        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName)

        ctrl = DTcontrol_C()

        # Run identifier
        ctrl.title = "test_ParticleMigration.py:test_2D_particle_migration"
        # Run author
        ctrl.author = "tph"
        
        ctrl.time = 0.0
        ctrl.timeloop_count = 0

        ctrl.dt = 0.5
        ctrl.n_timesteps = 19
        ctrl.MAX_FACET_CROSS_COUNT = 100

        # Run for 2D mesh

        # 1. Attach the particle mesh to particle_P
        # 2. Attach the C++ particle movers.
        # 3. Compute the cell-neighbors and facet-normals for the particle movers.
        self.particle_P.initialize_particle_mesh(self.pmesh2D)

        # Attach the particle integrators
        self.particle_P.initialize_particle_integration()

        ### Particle boundary-conditions

        # See the TODO of 4apr20 for redoing this.
        
        # Make a dictionary associating the above-named boundaries of the particle mesh with
        # user-supplied call-back functions.
        
        # First, we need to make a UserParticleBoundaryFunctions_... object, which
        # contains the boundary-crossing callback functions.
        
        spNames = self.particle_P.species_names
        
        # Import C++ particle boundary-conditions
        userParticleBoundaryFunctionsSOlibName = "user_particle_boundary_functions_cartesian_xyz_solib"
        if userParticleBoundaryFunctionsSOlibName not in sys.modules:        
            infoMsg = "%s\t\"Importing %s\"" % (fncName, userParticleBoundaryFunctionsSOlibName)
            print(infoMsg)
        userParticleBoundaryFunctionsSOlib = im_m.import_module(userParticleBoundaryFunctionsSOlibName)
        # Call the constructor to make a UserParticleBoundaryFunctions_... object
        userPBndFns = userParticleBoundaryFunctionsSOlib.UserParticleBoundaryFunctions(self.particle_P.position_coordinates)
        # Create the map from mesh facets to particle callback functions:
        pmeshBCs = self.particle_P.particle_solib.ParticleMeshBoundaryConditions(spNames, self.particle_P.pmesh_M, userPBndFns, print_flag=False)

        # Add pmeshBCs to the Particle_C object
        self.particle_P.pmesh_bcs = pmeshBCs
        
        # cell_dict{} is needed by the Python version of is_inside_cell(), which is
        # used in compute_mesh_cell_indices() below.
        self.particle_P.pmesh_M.compute_cell_dict()
        
        # Get the initial cell index of each particle.
        self.particle_P.compute_mesh_cell_indices()

        # Save the initial conditions (p_ic) for plotting
        p_ic = []
        sp = self.particle_P.neutral_species[0]
        for ip in [0, 1]:
            (pseg, offset) = self.particle_P.sap_dict[sp].get_segment_and_offset(ip)
            p = pseg[offset].copy()  # Have to make a copy! Otherwise you overwrite the only copy of the particle
            
#            p_ic.append(self.particle_P.sap_dict[sp].get(ip)) # Don't do this: it's a reference.
            p_ic.append(p)
#            print 'ip =', ip, 'p_ic =', p_ic[ip]

        # The expected final position and cell index

        # First particle

        xsp0 = -9.5; ysp0 =  -9.5; zsp0 = 0.0
        vxsp0 = -2.0; vysp0 = 0.0; vzsp0 = 0.0

        weight0 = 2.0
        bitflag0 = 2
        cell_index0 = 0
        unique_ID0 = 0
        crossings = 0

        psp0 = (xsp0,ysp0,zsp0, vxsp0,vysp0,vzsp0, weight0, bitflag0, cell_index0, unique_ID0, crossings)

        # Second particle

        xsp1 = -9.5; ysp1 =  -9.5; zsp1 = -9.5
        vxsp1 = -2.0; vysp1 = -2.0; vzsp1 = -2.0

        weight1 = 3.0
        bitflag1 = 2
        cell_index1 = 0
        unique_ID1 = 0
        crossings = 0

        psp1 = (xsp1,ysp1,zsp1, vxsp1,vysp1,vzsp1, weight1, bitflag1, cell_index1, unique_ID1, crossings)

        p_expected = (psp0, psp1)

        # Integrate for n_timesteps
        print("Advancing", self.particle_P.get_total_particle_count(), "particles for", ctrl.n_timesteps, "timesteps on a 2D mesh")
        for istep in range(ctrl.n_timesteps):
            self.particle_P.advance_neutral_particles(ctrl)

        # Create a mesh plotter to display the trajectory
        if self.plot_results is True:
          plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name + ": First & last positions"
          plotter=df_m.plot(self.particle_P.pmesh_M.mesh, title=plotTitle)

        # Check the results
        ncoords = self.particle_P.particle_dimension # number of particle coordinates to check
        for sp in self.particle_P.neutral_species:
            for ip in [0, 1]:
                (pseg, offset) = self.particle_P.sap_dict[sp].get_segment_and_offset(ip)
                getparticle = pseg[offset] # Retrieve the particle from the SAP.
                if self.plot_results is True:
                    #mplot_m.plot(data_arr['x'], data_arr['y'])
                    mplot_m.plot([p_ic[ip][0], getparticle[0]], [p_ic[ip][1], getparticle[1]])
#                print 'sp =', sp, 'expected =', p_expected[ip]
#                print 'calculated = ', getparticle
#                path = np_m.array([p_ic[ip][0], p_ic[ip][1], getparticle[0], getparticle[1]])
# Replace with point plot
#                plotter.add_polygon(path)

                for ic in range(ncoords):
                    self.assertAlmostEqual(p_expected[ip][ic], getparticle[ic], places=6, msg="Particle is not in correct position")
                cell_index_position = -3
#                print fncName, "expected cell =", p_expected[ip][cell_index_position], "computed cell =", getparticle[cell_index_position]
                self.assertEqual(p_expected[ip][cell_index_position], getparticle[cell_index_position], msg="Particle is not in correct cell")

        if self.plot_results is True:
            mplot_m.show()

        return
#    def test_2D_particle_migration(self):ENDDEF


    def test_3D_particle_migration(self):
        """Test particle migration across cells on a 3D mesh.

        """
        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName)

        ctrl = DTcontrol_C()

        # Run identifier
        ctrl.title = "test_ParticleMigration.py:test_3D_particle_migration"
        # Run author
        ctrl.author = "tph"
        
        ctrl.time = 0.0
        ctrl.timeloop_count = 0

        ctrl.dt = 0.5
        ctrl.n_timesteps = 19
        ctrl.MAX_FACET_CROSS_COUNT = 100

        # Run for 3D mesh

        # 1. Attach the particle mesh to particle_P
        # 2. Attach the C++ particle movers.
        # 3. Compute the cell-neighbors and facet-normals for the particle movers.
        self.particle_P.initialize_particle_mesh(self.pmesh3D)
        
        # Attach the particle integrators
        self.particle_P.initialize_particle_integration()

        ### Particle boundary-conditions

        # See the TODO of 4apr20 for redoing this.
        
        # Make a dictionary associating the above-named boundaries of the particle mesh with
        # user-supplied call-back functions.
        
        # First, we need to make a UserParticleBoundaryFunctions_... object, which
        # contains the boundary-crossing callback functions.
        
        spNames = self.particle_P.species_names
        
        # Import C++ particle boundary-conditions
        userParticleBoundaryFunctionsSOlibName = "user_particle_boundary_functions_cartesian_xyz_solib"
        if userParticleBoundaryFunctionsSOlibName not in sys.modules:
            infoMsg = "%s\t\"Importing %s\"" % (fncName, userParticleBoundaryFunctionsSOlibName)
            print(infoMsg)
        userParticleBoundaryFunctionsSOlib = im_m.import_module(userParticleBoundaryFunctionsSOlibName)
        # Call the constructor to make a UserParticleBoundaryFunctions_... object
        userPBndFns = userParticleBoundaryFunctionsSOlib.UserParticleBoundaryFunctions(self.particle_P.position_coordinates)
        # Create the map from mesh facets to particle callback functions:
        pmeshBCs = self.particle_P.particle_solib.ParticleMeshBoundaryConditions(spNames, self.particle_P.pmesh_M, userPBndFns, print_flag=False)        

        # Add pmeshBCs to the Particle_C object
        self.particle_P.pmesh_bcs = pmeshBCs

        # cell_dict{} is needed by the Python version of is_inside_cell(), which is
        # used in compute_mesh_cell_indices() below.
        self.particle_P.pmesh_M.compute_cell_dict()
        
        # Get the initial cell index of each particle.
        self.particle_P.compute_mesh_cell_indices()

        # Save the initial conditions (p_ic) for plotting
        p_ic = []
        sp = self.particle_P.neutral_species[0]
        for ip in [0, 1]:
            (pseg, offset) = self.particle_P.sap_dict[sp].get_segment_and_offset(ip)
            p = pseg[offset].copy()  # Have to make a copy! Otherwise you overwrite the only copy of the particle
            p_ic.append(p)

        # The expected final position and cell index

        # First particle

        xsp0 = -9.5; ysp0 =  -9.5; zsp0 = 0.0
        vxsp0 = -2.0; vysp0 = 0.0; vzsp0 = 0.0

        weight0 = 2.0
        bitflag0 = 2
        cell_index0 = 98
        unique_ID0 = 0
        crossings = 0

        psp0 = (xsp0,ysp0,zsp0, vxsp0,vysp0,vzsp0, weight0, bitflag0, cell_index0, unique_ID0, crossings)

        # Second particle

        xsp1 = -9.5; ysp1 =  -9.5; zsp1 = -9.5
        vxsp1 = -2.0; vysp1 = -2.0; vzsp1 = -2.0

        weight1 = 3.0
        bitflag1 = 2
        cell_index1 = 0
        unique_ID1 = 0
        crossings = 0

        psp1 = (xsp1,ysp1,zsp1, vxsp1,vysp1,vzsp1, weight1, bitflag1, cell_index1, unique_ID1, crossings)

        p_expected = (psp0, psp1)

        # Integrate for n_timesteps
        print("Advancing", self.particle_P.get_total_particle_count(), "particles for", ctrl.n_timesteps, "timesteps on a 3D mesh")
        for istep in range(ctrl.n_timesteps):
            self.particle_P.advance_neutral_particles(ctrl)

        # Create a mesh plotter to display the trajectory (just the
        # first and last positions)
        if self.plot_results is True:        
            plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name + ": First & last positions"
            plotter=df_m.plot(self.particle_P.pmesh_M.mesh, title=plotTitle)
        
        # Check the results
        ncoords = self.particle_P.particle_dimension # number of particle coordinates to check
        for sp in self.particle_P.neutral_species:
            for ip in [0, 1]:
                (pseg, offset) = self.particle_P.sap_dict[sp].get_segment_and_offset(ip)
                getparticle = pseg[offset] # Retrieve the particle from the SAP.
                if self.plot_results is True:
                    mplot_m.plot([p_ic[ip][0], getparticle[0]], [p_ic[ip][1], getparticle[1]], [p_ic[ip][2], getparticle[2]])
#                print 'expected = ', p_expected[ip]
#                print 'calculated = ', getparticle
#                path = np_m.array([p_ic[ip][0], p_ic[ip][1], p_ic[ip][2], getparticle[0], getparticle[1], getparticle[2]])
# Replace this with a point plot:
#                plotter.add_polygon(path)

                for ic in range(ncoords):
                    self.assertAlmostEqual(p_expected[ip][ic], getparticle[ic], places=6, msg="Particle is not in correct position")
                cell_index_position = -3
#                print (fncName, "expected cell =", p_expected[ip][cell_index_position], ", computed cell =", getparticle[cell_index_position])
                self.assertEqual(p_expected[ip][cell_index_position], getparticle[cell_index_position], msg="Particle is not in correct cell")

        if self.plot_results is True:
            mplot_m.show()

#        yesno = raw_input("Just called show() in test_3D_particle_migration")

        return
#    def test_3D_particle_migration(self):ENDDEF


#class TestCppParticleMigration:ENDCLASS

if __name__ == '__main__':
    unittest.main()
