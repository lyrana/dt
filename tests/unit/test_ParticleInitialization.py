#!/usr/bin/env python

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import numpy
import importlib as im_M
import unittest

from UserUnits_Module import MyPlasmaUnits_C
from Particle_Module import *

#STARTCLASS
class TestParticleInitialization(unittest.TestCase):
    """Test classes in Particle_Module"""
    
    def setUp(self):

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

        # Initializations performed before each test go here...

        # Create an instance of the DTparticleInput class
        pinCI = ParticleInput_C()
        # Set up particle variables
        pinCI.precision = numpy.float64

        pinCI.particle_integration_loop = 'loop-on-particles'
        pinCI.position_coordinates = ['x', 'y', 'z'] # determines the particle-storage dimensions
        pinCI.force_components = ['x', 'y','z']
        pinCI.force_precision = numpy.float64


        ### Particle species input

        # Give the properties of the particle species.  The charges
        # and masses are normally those of the physical particles, and
        # not the computational macroparticles.  Macroparticle weights
        # are specified or computed in a separate file (see
        # user_particle_module below) giving the distribution
        # functions, and can vary from particle to particle.

#         pinCI.particle_species = (('plasmaelectrons',
#                              {'initial_distribution_type' : 'listed',
#                               'charge' : -1.0*MyPlasmaUnits_C.elem_charge,
#                               'mass' : 1.0*MyPlasmaUnits_C.electron_mass,
#                               'dynamics' : 'implicit',
# #                              'number_per_cell' : 12,
#                               }
#                              ),
#                             ('Hplus', 
#                              {'initial_distribution_type' : 'listed',
#                               'charge' : 1.0*MyPlasmaUnits_C.elem_charge,
#                               'mass' : 1.0*MyPlasmaUnits_C.AMU,
#                               'dynamics' : 'implicit',
# #                              'number_per_cell' : 6,
#                               }
#                              ),
#                             ('He', 
#                              {'initial_distribution_type' : None,
#                               'charge' : 0.0,
#                               'mass' : 4.0*MyPlasmaUnits_C.AMU,
#                               'dynamics' : 'implicit',
# #                              'number_per_cell' : 1,
#                               }
#                              ),
#                             )

        speciesName = 'plasma_electrons'
        charge = -1.0*MyPlasmaUnits_C.elem_charge
        mass = 1.0*MyPlasmaUnits_C.electron_mass
        dynamics = 'implicit'
        plasmaElectronsCI = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        speciesName = 'H_plus'
        charge = 1.0*MyPlasmaUnits_C.elem_charge
        mass = 1.0*MyPlasmaUnits_C.AMU
        dynamics = 'implicit'
        HPlusCI = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        speciesName = 'He'
        charge = 0.0
        mass = 4.0*MyPlasmaUnits_C.AMU
        dynamics = 'implicit'
        HeCI = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        # Add these species to particle input
        pinCI.particle_species = (plasmaElectronsCI, HPlusCI, HeCI,
                                 )
        # Make the particle object from pinCI
        self.particleCI = Particle_C(pinCI, print_flag=False)

        # Give the name of the .py file containing special particle data (lists of
        # particles, boundary conditions, source regions, etc.)
        userParticleModule = "UserParticles_H_He_e"

        # Import this module
        UPrt_M = im_M.import_module(userParticleModule)

        self.particleCI.user_particle_module = userParticleModule
        self.particleCI.user_particle_class = userParticleClass = UPrt_M.UserParticleDistributions_C

        ### Provide input for particles present at t=0

        ## plasma electrons are present at t=0
        # Name the initialized species (it should be in species_names above)
        speciesName = 'plasma_electrons'
        # Check that this species has been defined above
        if speciesName not in self.particleCI.species_names:
            print fncName + "The species", speciesName, "has not been defined"
            sys.exit()

        # Specify how the species will be initialized
        initialDistributionType = 'listed'
        # Check that there's a function listing the particles particles
        printFlag = True
        if hasattr(userParticleClass, speciesName):
            if printFlag: print fncName + "(DnT INFO) Initial distribution for", speciesName, "is the function of that name in", userParticleClass
        # Write error message and exit if no distribution function exists
        else:
            errorMsg = fncName + "(DnT ERROR) Need to define a particle distribution function %s in %s for species %s " % (speciesName, userParticleModule, speciesName)
            sys.exit(errorMsg)

        # Collect the parameters into a dictionary
        # The 'listed' type will expect a function with the same name as the species.
        plasmaElectronParams = {'species_name': speciesName,
                              'initial_distribution_type': initialDistributionType,
                              }

        ## H+ ions are present at t=0
        speciesName = 'H_plus'
        # Check that this species has been defined above
        if speciesName not in self.particleCI.species_names:
            print fncName + "The species", speciesName, "has not been defined"
            sys.exit()

        # Specify how the species will be initialized
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
        HPlusParams = {'species_name': speciesName,
                              'initial_distribution_type': initialDistributionType,
                       }

        # The dictionary keys are mnemonics for the initialized particles
        initialParticlesDict = {'initial_plasma_electrons': (plasmaElectronParams,),
                                'initial_H_plus': (HPlusParams,),
                                }

        self.particleCI.initial_particles_dict = initialParticlesDict

        self.pinCI = pinCI # Used for the tests below

        return
#    def setUp(self):ENDDEF

#class TestParticleInitialization(unittest.TestCase):
    def test_species_names(self):
        """ Check that the species names are those specified by the user.
        """

        fncname = sys._getframe().f_code.co_name
        print '\ntest:', fncname, '('+__file__+')'

        particle_species = self.pinCI.particle_species
        # Check the names of the species
        for i in range(len(particle_species)):
            expected_name = particle_species[i].name
#            print "expected_name:", expected_name,"actual name:", self.particleCI.species_names[i]
            self.assertEqual(self.particleCI.species_names[i], expected_name, msg = "Species name is not correct")
        return
#    def test_species_names(self):ENDDEF

#class TestParticleInitialization(unittest.TestCase):
    def test_listed_particles(self):
        """ 1. Check that the stored particles have the values listed by the user.
            2. Check the methods that count the particles.
        """

        fncname = sys._getframe().f_code.co_name
        print '\ntest:', fncname, '('+__file__+')'

        user_particle_class = self.particleCI.user_particle_class

        ninput_total = 0

        initialParticlesDict = self.particleCI.initial_particles_dict

        # Loop on initialized particles
        for ip in initialParticlesDict:
            ipList = initialParticlesDict[ip]
            ipParams = ipList[0]
            s = ipParams['species_name']
            initialDistributionType = ipParams['initial_distribution_type']

            if initialDistributionType == 'listed':
                # Put user-listed particles into the storage array
                self.particleCI.create_from_list(s, False)
                # Get the original input data by calling the
                # user-provided function in 'user_particle_file'.
                # The function name is the species name, and the
                # argument for the particles used in this particular
                # test is 'listed'
                ninput, particles = getattr(user_particle_class, s)(initialDistributionType)
                ninput_total += ninput
                # Check that the particles has the user-input values
                for i in range(ninput):
                    putparticle = particles[i]
                    # Get the stored particle data
#                    getparticle = self.particles.pseg_arr[sp][ip]
                    getparticle = self.particleCI.pseg_arr[s].get(i)
#                    print 'putparticle = ', putparticle
#                    print 'getparticle = ', getparticle
                    for ix in range(len(getparticle)):
                        self.assertAlmostEqual(getparticle[ix], putparticle[ix], msg="Particle is not correct")

                # Check the number of stored particles in each species
                nstored = self.particleCI.get_species_particle_count(s, print_flag = False)
                self.assertEqual(nstored, ninput, msg = "Number of stored particles is not correct")

        # check total number of stored particles
        nstored = self.particleCI.get_total_particle_count(print_flag = False)
        self.assertEqual(nstored, ninput_total, msg = "Number of stored particles is not correct")

        return
#    def test_listed_particles(self):ENDDEF

#class TestParticleInitialization(unittest.TestCase):
    def test_something(self):
        """ Check the number of particles in each species.
        """
        fncname = sys._getframe().f_code.co_name
        print '\ntest:', fncname, '('+__file__+')'

        return
#    def test_something(self):ENDDEF

if __name__ == '__main__':
    unittest.main()
