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

# The C++ functions are bound to Python names in a .so
#import dolfin_cpp
#import p_cpp_cartesian_xyz

#STARTCLASS
class TestCppParticleMigration(unittest.TestCase):
    """Test use of pybind11 to interface with C++

       Use a C++ function to advance particles, passing needed parameters from
       Python.  Neutral helium is the only species.  Initial particles are specified
       in UserParticles_3D.py.

    """
    
    def setUp(self):

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

        # Initializations performed before each test go here...

        # Create an instance of the DTparticleInput class
        pin = ParticleInput_C()

        # Set up particle variables
        pin.precision = np_m.float64
        pin.particle_integration_loop = 'loop-on-particles'

        # The particle coordinate system is 3D Cartesian
        pin.coordinate_system = 'cartesian_xyz'

        # Use the 3 position coordinates, since we're doing 1, 2, and 3D particle motion
        # This could be derived from the coordinate system:

# Moved this to Particle_C.init, based on coordinate system        
#        pin.position_coordinates = ['x', 'y', 'z'] # Determines particle storage dimension

        
        # Neutral particles: No forces.
        """
        pin.force_components = ['x', 'y',]
        """
        pin.force_precision = np_m.float64
        pin.use_cpp_movers = True        

        # Give the properties of the particle species.  The charges and masses are
        # normally those of the physical particles, and not the computational
        # macroparticles.  Macroparticle weights are specified or computed in a
        # separate file (see user_particles_module_name below) specifying the
        # particle distribution functions, and can vary from particle to particle.

        speciesName = 'neutral_H'
        charge = 0.0
        mass = 1.0*MyPlasmaUnits_C.AMU
        dynamics = 'neutral'
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
        self.pmesh1D = UserMesh_C(umi1D, compute_dictionaries=True, compute_tree=True, plot_flag=plotFlag, plot_title=plotTitle + ": 1D")

# Individual geometry dictionaries:
#        self.pmesh1D.compute_cell_vertices_dict()
#        self.pmesh1D.compute_cell_dict()


        return
#    def setUp(self):ENDDEF

    def test_1D_particle_migration(self):
        """Test that we can call a C++ function to push neutral particles.

           Initial particle positions are in UserParticles_3D.py.

           The first particle moves in -x only, from 9.5 to -9.5.
           The first particle starts at:
               (x0, y0, z0) = (9.5, -9.5, 0.0), with velocity:
               (ux0, uy0, uz0) = (-2.0, 0.0, 0.0)
           It moves in -x only, to -9.5

           The second particle starts at
               (x1, y1, z1) = (9.5, 9.5, 9.5)
               (ux1, uy1, uz1) = (-2.0, -2.0, -2.0)
           It moves in -x, -y, -z to the opposite corner.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName)

        ctrl = DTcontrol_C()

        # Run identifier
        ctrl.title = "ParticleAdvance using C++"
        # Run author
        ctrl.author = "tph"

        ctrl.time = 0.0
        ctrl.timeloop_count = 0

        ctrl.dt = 0.5
        ctrl.n_timesteps = 19
        ctrl.MAX_FACET_CROSS_COUNT = 100

        # Run on a 1D mesh
        # Note: this would be a copy of the field mesh for parallel runs.
        #self.particle_P.pmesh_M = self.pmesh1D
        # Need to compute the mesh facet-normals for particle movers
        
        self.particle_P.initialize_particle_mesh(self.pmesh1D);
        
#        self.particle_P.pmesh_M.compute_cell_facet_normals_dict(use_cpp=True)
#    def compute_cell_facet_normals_dict(self, use_cpp=False):

        # Get the initial cell index of each particle.
        # Note: this could be moved down closer to the beginning of the particle-advance loop.
        self.particle_P.compute_mesh_cell_indices()

        ### Put the expected ending results for the two particles into the p_expected tuple ###

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
        unique_ID1 = 1
        crossings = 0

        psp1 = (xsp1,ysp1,zsp1, vxsp1,vysp1,vzsp1, weight1, bitflag1, cell_index1, unique_ID1, crossings)

        p_expected = (psp0, psp1)

        #
        # Move the particles and check the final positions
        #
        
        speciesName = 'neutral_H'
        psa = self.particle_P.pseg_arr[speciesName] # segmented array for this species

        # Print the cell indices
#        print("cell_vertices_dict = ", self.particle_P.pmesh_M.cell_vertices_dict)

        # Integrate for n_timesteps
        print("Moving", self.particle_P.get_total_particle_count(), "particles for", ctrl.n_timesteps, "timesteps")
        for istep in range(ctrl.n_timesteps):
#            self.particle_P.move_neutral_particles_CPP(ctrl)
            print(fncName, "istep:", istep)
            self.particle_P.move_neutral_particles(ctrl)

        # Check the results
        ncoords = self.particle_P.particle_dimension # number of particle coordinates to check
        for sp in self.particle_P.neutral_species:
            for ip in [0, 1]:
                # getparticle = self.particle_P.pseg_arr[sp].get(ip)
                # Instead of .get(), retrieve the particle structure using the returned Numpy array it's in.
                (pseg, offset) = self.particle_P.pseg_arr[sp].get_segment_and_offset(ip)
                getparticle = pseg[offset] # Retrieve the particle from the SAP.
                print('expected = ', p_expected[ip])
                print('calculated = ', getparticle)
                for ic in range(ncoords):
                    self.assertAlmostEqual(p_expected[ip][ic], getparticle[ic], places=6, msg="Particle is not in correct position")
                cell_index_position = -3
#                print fncName, "expected cell =", p_expected[ip][cell_index_position], "computed cell =", getparticle[cell_index_position]
                self.assertEqual(p_expected[ip][cell_index_position], getparticle[cell_index_position], msg="Particle is not in correct cell")

        return
    # ENDDEF: def test_1D_particle_migration(self)


#class TestCppParticleMigration:ENDCLASS

if __name__ == '__main__':
    unittest.main()
