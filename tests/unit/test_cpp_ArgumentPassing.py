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

import test_cpp

#STARTCLASS
class TestPybind11(unittest.TestCase):
    """Test use of pybind11 to interface with C++.

       Test that we can pass and access attributes of a DnTcontrol_C object in C++.
       Test that we can pass and access attributes of Particle_C object in C++.

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
        pin.coordinate_system = 'cartesian_x_y_z'

        # Use the 3 position coordinates, since we're doing 1, 2, and 3D particle motion
        # This could be derived from the coordinate system:

# Moved this to Particle_C.init, based on coordinate system        
#        pin.position_coordinates = ['x', 'y', 'z'] # Determines particle storage dimension

        
        # Neutral particles: No forces.
        """
        pin.force_components = ['x', 'y',]
        """
        pin.force_precision = np_m.float64

        # Give the properties of the particle species.  The charges and masses are
        # normally those of the physical particles, and not the computational
        # macroparticles.  Macroparticle weights are specified or computed in a
        # separate file (see user_particles_module_name below) specifying the
        # particle distribution functions, and can vary from particle to particle.

        speciesName = 'neutral_H'
        charge = 1.0
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
        # Particles have a 3D/3D phase-space even though mesh is just 2D
        userParticlesModuleName = "UserParticles_3D"
#        print "setUp: UserParticle module is", userParticlesModuleName

        # Import this module
        userParticlesModule = im_m.import_module(userParticlesModuleName)

        # self.particle_P.user_particles_module_name = userParticlesModuleName
        self.particle_P.user_particles_class = userParticlesClass = userParticlesModule.UserParticleDistributions_C

        ### neutral H atoms are present at t=0
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

        # Create the initial particles
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

        # Create 1D, 2D and 3D meshes that the particles can be tested against

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
#        self.pmesh1D.compute_cell_vertex_dict()
#        self.pmesh1D.compute_cell_dict()

        return
#    def setUp(self):ENDDEF

    def test_1_pass_DnTcontrol(self):
        """Test that we can pass and access attributes of a DnTcontrol Python object
           in C++.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName)

        ctrl = DTcontrol_C()

        # Run identifier
        ctrl.title = "Test of argument-passing from Python to C++"
        # Run author
        ctrl.author = "tph"

        ctrl.time = 0.0
        ctrl.dt = 0.5
        ctrl.n_timesteps = 19

        # Pass a DT_control argument to C++:
        test_cpp.function_with_DTcontrol_arg(ctrl)

        return
#    def test_1_pass_DnTcontrol:ENDDEF

    def test_2_pass_mesh_M(self):
        """Test that we can pass and access attributes of Particle_C object.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName)

        # Run on a 1D mesh
        self.particle_P.pmesh_M = self.pmesh1D

        # Get the initial cell index of each particle.
        self.particle_P.compute_mesh_cell_indices()

        ### Put the expected ending results into the p_expected tuple ###

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

# Test passing the 1D particle storage to C++

        speciesName = 'neutral_H'
        psa = self.particle_P.pseg_arr[speciesName] # segmented array for this species
        
        test_cpp.function_with_particle_P_arg(self.particle_P)

        return
#    def test_2_pass_mesh_M:ENDDEF

#class TestPybind11:ENDCLASS

if __name__ == '__main__':
    unittest.main()
