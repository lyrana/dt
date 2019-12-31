#!/usr/bin/env python3

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import os
import importlib as im_m
import unittest

import numpy as np_m

import dolfin as df_m

from DT_Module import DTcontrol_C

from Dolfin_Module import *
from Particle_Module import *
from UserMesh_y_Fields_FE_XYZ_Module import *

from UserUnits_Module import MyPlasmaUnits_C

#STARTCLASS
class TestParticleInitialization(unittest.TestCase):
    """Test classes in Particle_Module"""
    
    def setUp(self):

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

        # Initializations performed before each test go here...

        # Create an instance of the DTparticleInput class
        pin = self.pin = ParticleInput_C()
        # Set up particle variables
        pin.precision = np_m.float64

        pin.particle_integration_loop = 'loop-on-particles'
        pin.coordinate_system = 'cartesian_xyz'
        pin.force_components = ['x', 'y', 'z']
        pin.force_precision = np_m.float64


        ### Particle species input

        # Give the properties of the particle species.  The charges
        # and masses are normally those of the physical particles, and
        # not the computational macroparticles.  Macroparticle weights
        # are specified or computed in a separate file (see
        # user_particles_module_name below) giving the distribution
        # functions, and can vary from particle to particle.

        speciesName = 'plasma_electrons'
        charge = -1.0*MyPlasmaUnits_C.elem_charge
        mass = 1.0*MyPlasmaUnits_C.electron_mass
        dynamics = 'implicit'
        plasmaElectrons_PS = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        speciesName = 'H_plus'
        charge = 1.0*MyPlasmaUnits_C.elem_charge
        mass = 1.0*MyPlasmaUnits_C.AMU
        dynamics = 'implicit'
        Hplus_PS = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        speciesName = 'He'
        charge = 0.0
        mass = 4.0*MyPlasmaUnits_C.AMU
        dynamics = 'implicit'
        He_PS = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        # Add these 3 species to particle input
        pin.particle_species = (plasmaElectrons_PS, Hplus_PS, He_PS,
                                 )
        # Make the particle object from pin
        p_P = self.particle_P = Particle_C(pin, print_flag=False)

        # Give the name of the .py file containing special particle data (lists
        # of particles, boundary conditions, particle-initialization regions,
        # particle source regions, etc.)
        userParticlesModuleName = "UserParticles_H_He_e"

        # Import this module
        userParticlesModule = im_m.import_module(userParticlesModuleName)

        # p_P.user_particles_module_name = userParticlesModuleName
        p_P.user_particles_class = userParticlesClass = userParticlesModule.UserParticleDistributions_C

        ### Provide input for particles present at t=0

        ## plasma electrons are present at t=0
        # Name the initialized species (it should be in species_names above)
        speciesName = 'plasma_electrons'
        # Check that this species has been defined above
        if speciesName not in p_P.species_names:
            print(fncName + "The species", speciesName, "has not been defined")
            sys.exit()

        # Specify how the species will be initialized
### Could you get this info without using this variable?
### E.g., by giving the name of the function that lists the particles?
        initialDistributionType = 'listed'
        # Check that there's a function listing the particles particles
        printFlag = True
        if hasattr(userParticlesClass, speciesName):
            if printFlag: print(fncName + "(DnT INFO) Initial distribution for", speciesName, "is the function of that name in", userParticlesClass)
        # Write error message and exit if no distribution function exists
        else:
            errorMsg = fncName + "(DnT ERROR) Need to define a particle distribution function %s in %s for species %s " % (speciesName, userParticlesModuleName, speciesName)
            sys.exit(errorMsg)

        # Collect the parameters into a dictionary
        # The 'listed' type will expect a function with the same name as the species.

### Here, we'd provide the name of the distribution function, 'plasma_electrons_listed'

        plasmaElectronParams = {'species_name': speciesName,
                                'initial_distribution_type': initialDistributionType,
                              }

        ## H+ ions are present at t=0
        speciesName = 'H_plus'
        # Check that this species has been defined above
        if speciesName not in p_P.species_names:
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
        HplusParams = {'species_name': speciesName,
                       'initial_distribution_type': initialDistributionType,
                       }

        # The dictionary keys are mnemonics for the initialized particles
        initialParticlesDict = {'initial_plasma_electrons': (plasmaElectronParams,),
                                'initial_H_plus': (HplusParams,),
                                }

        p_P.initial_particles_dict = initialParticlesDict

        return
#    def setUp(self):ENDDEF

#class TestParticleInitialization(unittest.TestCase):
    def test_1_species_names(self):
        """ Check that the species names are those specified by the user.
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest:', fncName, '('+__file__+')')

        particle_species = self.pin.particle_species
        # Check the names of the species
        for i in range(len(particle_species)):
            expected_name = particle_species[i].name
#            print "expected_name:", expected_name,"actual name:", self.particle_P.species_names[i]
            self.assertEqual(self.particle_P.species_names[i], expected_name, msg = "Species name is not correct")
        return
#    def test_species_names(self):ENDDEF

#class TestParticleInitialization(unittest.TestCase):
    def test_2_listed_particles(self):
        """ 1. Check that the stored particles have the values listed by the user.
            2. Check the methods that count the particles.
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName, '('+__file__+')')

        userParticlesClass = self.particle_P.user_particles_class

        ninput_total = 0

        initialParticlesDict = self.particle_P.initial_particles_dict

        # Loop on initialized particles
        for ip in initialParticlesDict:
            ipList = initialParticlesDict[ip]
            ipParams = ipList[0]
            s = ipParams['species_name']
            initialDistributionType = ipParams['initial_distribution_type']

            if initialDistributionType == 'listed':
                # Put user-listed particles into the storage array
                self.particle_P.create_from_list(s, False)
                # Get the original listed-particle input data from the
                # user-provided function in 'user_particles_file'.
                # The function name is the species name, and the
                # argument for the particles used in this particular
                # test is 'listed'
                ninput, particles = getattr(userParticlesClass, s)(initialDistributionType)
                ninput_total += ninput
                # Check that the particles stored in the SA have the user-input values
                for i in range(ninput):
                    putparticle = particles[i]
                    # Get the stored particle data
                    (pseg, offset) = self.particle_P.pseg_arr[s].get_segment_and_offset(i)
                    getparticle = pseg[offset] # Retrieve the particle from the SAP.
                    
#                    print 'putparticle = ', putparticle
#                    print 'getparticle = ', getparticle
                    for ix in range(len(getparticle)):
                        self.assertAlmostEqual(getparticle[ix], putparticle[ix], msg="Particle is not correct")

                # Check the number of stored particles in each species
                nStored = self.particle_P.get_species_particle_count(s, print_flag = False)
                self.assertEqual(nStored, ninput, msg = "Number of stored particles is not correct")

        # check total number of stored particles
        nStored = self.particle_P.get_total_particle_count(print_flag = False)
        self.assertEqual(nStored, ninput_total, msg = "Number of stored particles is not correct")

        return
#    def test_listed_particles(self):ENDDEF

#class TestParticleInitialization(unittest.TestCase):
    def test_3_function_initialized_particles(self):
        """ Test initialization of particles using a function over a region of
            the mesh.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName, '('+__file__+')')

        ## Control struct
        ctrl = DTcontrol_C()

        ## Mesh input
        mi2d_UMI = UserMeshInput_C()

# Replace with TUPLES since we're at toplevel
        mi2d_UMI.pmin = df_m.Point(-10.0, -10.0)
        mi2d_UMI.pmax = df_m.Point(10.0, 10.0)
        mi2d_UMI.cells_on_side = (2, 2)
        mi2d_UMI.diagonal = 'left'

        ## Create a 2D particle mesh, since we're going to initialize particles
        ## on a region of a mesh.
        # UserMesh_y_Fields_FE_XYZ_Module can make the mesh from the above input.
        plotFlag = False
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name + ": mesh"
        pmesh2d_UM = UserMesh_C(mi2d_UMI, compute_dictionaries=True, compute_cpp_arrays=False, compute_tree=True, plot_flag=plotFlag, plot_title=plotTitle)

        p_P = self.particle_P

        userParticlesClass = p_P.user_particles_class

        ## 1. Initial hot electrons on bottom-left side of the mesh

        # a. Provide the species name and physical parameters
        #    for this particle-initialization
        speciesName = 'plasma_electrons'

        # Check that this species has been defined
        if speciesName in p_P.species_names:
            charge = p_P.charge[speciesName]
            mass = p_P.mass[speciesName]
        else:
            print("The species", speciesName, "needs to be defined.")
            sys.exit()

        initialDistributionType = 'function_over_region'
        # Specify the initial physical number-density for this species
        numberDensity = 1.0e12
        # Check for positivity
        if numberDensity <= 0:
            print("Check number density for species", speciesName, "is zero or negative. Should be positive")
            sys.exit()

        # Compute a value for thermalSpeed from the temperature in eV
        temperature_eV = 2.0
        temp_joule = temperature_eV*MyPlasmaUnits_C.elem_charge
        thermalSpeed = np_m.sqrt(2.0*temp_joule/mass)

        # Set a drift velocity
        # Set a drift velocity
        vdrift_1 = 1.0
        vdrift_2 = 2.0
        vdrift_3 = 3.0

        velocity_coordinate_system = 'x_y_z' # This sets the interpretation of the
                                             # following values
        driftVelocity = (vdrift_1, vdrift_2, vdrift_3)
#        driftVelocity[:] = (1.0,)

        # Set desired particles-per-cell
        numberPerCell = 1000

        # Collect the parameters for this initialization in a dictionary
        hotElectronParams = {'species_name': speciesName,
                             'initial_distribution_type': initialDistributionType,
                             'number_density': numberDensity,
                             'thermal_speed': thermalSpeed,
                             'velocity_coordinate_system': velocity_coordinate_system,
                             'drift_velocity': driftVelocity,
                             'number_per_cell': numberPerCell}

        # b. Name the particle-creation function.
        maxwellianGenerator = p_P.create_maxwellian_particles

        # c. Initialization-region geometry
        xmin = -9.0; ymin = -9.0
        xmax = -4.0; ymax = -4.0

# Don't use df_m in the toplevel routine, if possible
#        pmin = df_m.Point(xmin, ymin)
#        pmax = df_m.Point(xmax, ymax)
        pmin = (xmin, ymin)
        pmax = (xmax, ymax)
        hotElectronsRegion_RR = RectangularRegion_C(pmesh2d_UM, pmin, pmax)

        ## Put all the particle-initialization data into a dictionary

        # The dictionary keys are mnemonics for particle-initialization names.
        # The dictionary values are lists of tuples giving the 
        #   a. the physical initialization parameters.
        #   b. particle creation function and
        #   c. initialization region

        ## Note: the names here are particle-initialization names, NOT species names!
        initialParticlesDict = {'initial_hot_electrons': (hotElectronParams, maxwellianGenerator, hotElectronsRegion_RR),
#                              'initial_background_electrons': (backgroundElectronParams, maxwellianGenerator, wholeMesh),
                              }

        # Replace the 'listed' initial particles made in Setup() with
        # the 'function_initialized' initial particles created above

        p_P.initial_particles_dict = initialParticlesDict

        # Loop on the initialization methods in the dictionary

        # For 'listed' initialization, the dictionary has entries like
        # {ipName: (ipParams,)}
        # For 'function_over_region' initialization, the dictionary has entries
        # like
        # {ipName: (ipParams, ipFunc, ipRegion_RR)}

        nCreatedTotal = 0
        for ipName, ipTuple in initialParticlesDict.items():
            ipParams = ipTuple[0]
            s = ipParams['species_name']
            initialDistributionType = ipParams['initial_distribution_type']

            if initialDistributionType == 'function_over_region':
                print(fncName, "Initializating", ipName, "particles")
                ipFunc = ipTuple[1]
                ipRegion_RR = ipTuple[2]
                # Invoke the creation function
                step = 0
                time = 0.0
                neg_E_field = None
                ipFunc(step, time, ipRegion_RR, ipParams, neg_E_field)

                # Check the number of stored particles for each species

                nCreated = ipRegion_RR.ncell * ipParams['number_per_cell']
                print(fncName, "Number of particles created =", nCreated)
                nCreatedTotal += nCreated

                nStored = p_P.get_species_particle_count(s, print_flag = False)
                self.assertEqual(nStored, nCreated, msg = "Number of stored particles is not correct")

        # check total number of stored particles
        nStored = p_P.get_total_particle_count(print_flag = False)
        print(fncName, "Total number of particles created =", nCreatedTotal)
        self.assertEqual(nStored, nCreatedTotal, msg = "Total number of stored particles is not correct")

        # Write out the particles to an HDF5 file

        ctrl.timeloop_count = 0
        ctrl.time = 0.0

        # Run identifier
        ctrl.title = "test_ParticleGeneration"
        # Run author
        ctrl.author = "tph"

        ### Select output for particles
        ctrl.particle_output_file = "particleInitialization.h5part"
        ctrl.particle_output_interval = 1
        ctrl.particle_output_attributes = ('species_index', 'x', 'y', 'z', 'ux', 'uy', 'uz', 'weight')

        # Check these values
        p_P.check_particle_output_parameters(ctrl)

        # Dump the particle data to a file
        p_P.initialize_particle_output_file(ctrl)

        # Write the particle attributes
        p_P.write_particles_to_file(ctrl)

        return
#    def test_function_initialized_particles(self):ENDDEF


#class TestParticleInitialization(unittest.TestCase):
    def test_something(self):
        """ Check the number of particles in each species.
        """
        fncName = sys._getframe().f_code.co_name
        print('\ntest:', fncName, '('+__file__+')')

        return
#    def test_something(self):ENDDEF

if __name__ == '__main__':
    unittest.main()
