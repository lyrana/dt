#!/usr/bin/env python3

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import numpy
import importlib as im_m
import unittest

import DT_Module as DT_m

from UserUnits_Module import MyPlasmaUnits_C
from Particle_Module import *

import particle_cartesian_xyz_solib as p_cpp

#STARTCLASS
class Vec_C(object):
    """ Creates a 1, 2, or 3D vector.
    """

    def __init__(self, dimension, values):
        if dimension >= 1:
            self.x = values[0]
        if dimension >= 2:
            self.y = values[1]
        if dimension >= 3:
            self.z = values[2]

        return

# class Vec_C(object): ENDCLASS

#STARTCLASS
class TestParticleUniformE(unittest.TestCase):
    """Test classes in Particle_Module"""
    
    def setUp(self):

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

        # initializations for each test go here...

        # Create an instance of the ParticleInput class
        pin = ParticleInput_C()
        # Initialize particles
        pin.precision = numpy.float64
        pin.particle_integration_loop = 'loop-on-particles'
        pin.coordinate_system = 'cartesian_xyz'        
#        pin.position_coordinates = ['x', 'y', 'z'] # determines the particle-storage dimensions
        pin.force_components = ['x', 'y',]
        pin.force_precision = numpy.float64
        pin.use_cpp_movers = True
        
        # Specify the particle species for this calculation

        speciesName = 'plasma_electrons'
        charge = -1.0*MyPlasmaUnits_C.elem_charge
        mass = 1.0*MyPlasmaUnits_C.electron_mass
        dynamics = 'explicit'
        plasmaElectrons_S = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        speciesName = 'H_plus'
        charge = 1.0*MyPlasmaUnits_C.elem_charge
        mass = 1.0*MyPlasmaUnits_C.proton_mass
        dynamics = 'explicit'
        HPlus_S = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        speciesName = 'He'
        charge = 0.0
        mass = 4.0*MyPlasmaUnits_C.AMU
        dynamics = 'explicit'
        He_S = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        # Add these species to particle input
        pin.particle_species = (plasmaElectrons_S, HPlus_S, He_S,
                                )
        # Make the particle object from pin
        self.particle_P = Particle_C(pin, print_flag=False)

        # Give the name of the .py file containing additional particle data (lists of
        # particles, boundary conditions, source regions, etc.)
        userParticlesModuleName = "UserParticles_H_He_e"

        # Import this module
        userParticlesModule = im_m.import_module(userParticlesModuleName)

        #self.particle_P.user_particles_module_name = userParticlesModuleName
        self.particle_P.user_particles_class = userParticlesClass = userParticlesModule.UserParticleDistributions_C


        ### plasma_electrons are present at t=0

        # Name the initialized species (it should be in species_names above)
        speciesName = 'plasma_electrons'
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
            errorMsg = fncName + "(DnT ERROR) Need to define a particle distribution function %s in %s for species %s " % (speciesName, userParticlesModuleName, speciesName)
            sys.exit(errorMsg)

        # Collect the parameters into a dictionary
        # The 'listed' type will expect a function with the same name as the species.
        plasmaElectronParams = {'species_name': speciesName,
                              'initial_distribution_type': initialDistributionType,
                              }

        ## H+ ions are present at t=0
        speciesName = 'H_plus'
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
        HPlusParams = {'species_name': speciesName,
                              'initial_distribution_type': initialDistributionType,
                       }

        # The dictionary keys are mnemonics for the initialized particles
        initialParticlesDict = {'initial_plasma_electrons': (plasmaElectronParams,),
                                'initial_H_plus': (HPlusParams,),
                                }

        self.particle_P.initial_particles_dict = initialParticlesDict

        # Create the initial particles
        for ip in initialParticlesDict:
            ipList = initialParticlesDict[ip]
            ipParams = ipList[0]
            s = ipParams['species_name']
            initialDistributionType = ipParams['initial_distribution_type']
            if initialDistributionType == 'listed':
                # Put user-listed particles into the storage array
                self.particle_P.create_from_list(s, False)

        return
#    def setUp(self):ENDDEF

#class TestParticleUniformE(unittest.TestCase):
    def test_1_electric_field_push_1step(self):
        """ Check that the electric field push is correct.  Push
            sample particles for 1 step.
        """
        fncName = sys._getframe().f_code.co_name
        print('\ntest: ', fncName, '('+__file__+')')

        ctrl = DT_m.DTcontrol_C()

        ctrl.dt = 1.0e-5
        ctrl.n_timesteps = 1

        # Set constant fields
        E0 = (1.0e-4, 2.0e-4, 3.0e-4)
#        B0 = (0.0, 0.0, 0.0)
        ctrl.E0 = Vec_C(self.particle_P.particle_dimension, E0)
#        ctrl.B0 = Vec_C(self.particle_P.particle_dimension, B0)

        # The EXPECTED results, from ParticleUniformE.ods, one for each species:
        		
        # First species
        xsp1 = 0.0282411799; ysp1 = 0.0264823597; zsp1 = 0.0247235396
        vxsp1 = 2824.117985104; vysp1 = 1648.2359702079; vzsp1 = 472.3539553119
        weight1=1.0
        bitflag1 = 1
        psp1 = (xsp1,ysp1,zsp1, vxsp1,vysp1,vzsp1, weight1, bitflag1)

        # Second species
        xsp2 = 3.9578833918911e-06; ysp2 = 1.39157667837822e-05; zsp2 = 2.38736501756733e-05
        vxsp2 = 0.3957883392; vysp2 = 0.3915766784; vzsp2 = 0.3873650176
        weight2 = 3.0
        bitflag2 = 1
        psp2 = (xsp2,ysp2,zsp2, vxsp2,vysp2,vzsp2, weight2, bitflag2)
        p_expected = (psp1, psp2)

        # Loop on the species, moving the particles one timestep, and
        # then checking the particle end positions
        ncoords = self.particle_P.particle_dimension # number of particle coordinates to check
        isp = 0

#        print("Moving", self.particle_P.get_total_particle_count(), "particles for", ctrl.n_timesteps, "timesteps")
        for sp in self.particle_P.species_names:
            if self.particle_P.get_species_particle_count(sp) == 0: continue

            p_cpp.move_charged_species_in_uniform_fields(sp, ctrl, self.particle_P)
            
            # Check that the first particles in the array reach the right speed
#            getparticle = self.particle_P.pseg_arr[sp][0]
#            getparticle = self.particle_P.pseg_arr[sp].get(0)
            (pseg, offset) = self.particle_P.pseg_arr[sp].get_segment_and_offset(0)
            getparticle = pseg[offset]
            
#            print 'calculated = ', getparticle
#            print 'expected = ', p_expected[isp]
            for ic in range(ncoords):
#            for ix in range(len(getparticle)):
                self.assertAlmostEqual(getparticle[ic], p_expected[isp][ic], msg="Particle is not in correct position after 1 step")
            isp += 1

        return
#    def test_1_electric_field_push_1step(self):ENDDEF

#class TestParticleUniformE(unittest.TestCase):
    def test_2_electric_field_push_10steps(self):
        """ Check that the electric field push is correct.  Push
            sample particles for 10 steps.
        """
        fncName = sys._getframe().f_code.co_name
        print('\ntest: ', fncName, '('+__file__+')')

        ctrl = DT_m.DTcontrol_C()

        ctrl.dt = 1.0e-5
        ctrl.n_timesteps = 10

        # Set constant fields
        E0 = (1.0e-4, 2.0e-4, 3.0e-4)
#        B0 = (0.0, 0.0, 0.0)
        ctrl.E0 = Vec_C(self.particle_P.particle_dimension, E0)
#        ctrl.B0 = Vec_C(self.particle_P.particle_dimension, B0)

        # The expected results, from test-results.ods, one for each species:
        		
        # First species

        xsp1 = 0.2032648918; ysp1 =  0.0165297836; zsp1 = -0.1702053246
        vxsp1 = 1241.1798510396; vysp1 = -1517.6402979208; vzsp1 = -4276.4604468811

        weight1=1.0
        bitflag1 = 1

        psp1 = (xsp1,ysp1,zsp1, vxsp1,vysp1,vzsp1, weight1, bitflag1)

        # Second species

        xsp2 = 8.26835865540104E-005; ysp2 = 0.0001353672; zsp2 = 0.0001880508
        vxsp2 = 1.2578833919; vysp2 = 2.1157667838; vzsp2 = 2.9736501757

        weight2 = 3.0
        bitflag2 = 1

        psp2 = (xsp2,ysp2,zsp2, vxsp2,vysp2,vzsp2, weight2, bitflag2)

        p_expected = (psp1, psp2)

        # Integrate for n_timesteps
        ctrl.time_step = 0
        ctrl.time = 0.0
#        print("Moving", self.particle_P.get_total_particle_count(), "particles for", ctrl.n_timesteps, "timesteps")        
        for istep in range(ctrl.n_timesteps):
#            print 'test_2_electric_field_push_10steps: istep =', istep

            # Loop on the species
            for sp in self.particle_P.species_names:
#                if sp != 'plasmaelectrons': continue
#                print 'test_2_electric_field_push_10steps: sp =', sp
                if self.particle_P.get_species_particle_count(sp) == 0: continue

                p_cpp.move_charged_species_in_uniform_fields(sp, ctrl, self.particle_P)

            ctrl.time_step += 1
            ctrl.time += ctrl.dt

        # check the final position of the first particle of each species
        ncoords = self.particle_P.particle_dimension # number of particle coordinates to check
        isp = 0
        for sp in self.particle_P.species_names:
            if self.particle_P.get_species_particle_count(sp) == 0: continue
            # Check that the first particles in the array reach the right speed
#            getparticle = self.particle_P.pseg_arr[sp].get(0)
            (pseg, offset) = self.particle_P.pseg_arr[sp].get_segment_and_offset(0)
            getparticle = pseg[offset]
#            print('calculated = ', getparticle, 'for species', sp)
#            print('expected = ', p_expected[isp], 'for species', sp)
            for ic in range(ncoords):
#            for ix in range(len(getparticle)):
                self.assertAlmostEqual(getparticle[ic], p_expected[isp][ic], msg="Particle is not in correct position after 10 steps")
            isp += 1

        return
#    def test_2_electric_field_push_10steps(self):ENDDEF

#class TestParticleUniformE(unittest.TestCase):
    def test_3_magnetic_field_push(self):            
        """ 1. Check that the magnetic field push is correct.
            2. 
        """
        fncName = sys._getframe().f_code.co_name
        print('\ntest: ', fncName, '('+__file__+')')
        pass

        return

#class TestParticleUniformE(unittest.TestCase):
    def test_4_something(self):
        """ Check the number of particles in each species.
        """
        fncName = sys._getframe().f_code.co_name
        print('\ntest: ', fncName, '('+__file__+')')
        pass

        return

# class TestParticleUniformE(unittest.TestCase): ENDCLASS


if __name__ == '__main__':
    unittest.main()
