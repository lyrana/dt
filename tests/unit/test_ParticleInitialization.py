#!/usr/bin/env python

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import numpy
import importlib as im_M
import unittest

from DT_Module import DTparticleInput_C
from UserUnits_Module import MyPlasmaUnits_C
from Particle_Module import Particle_C

class TestParticleInitialization(unittest.TestCase):
    """Test classes in Particle_Module"""
    
    def setUp(self):

        # Initializations performed before each test go here...

        # Create an instance of the DTparticleInput class
        pinCI = DTparticleInput_C()
        # Set up particle variables
        pinCI.precision = numpy.float64

        pinCI.particle_integration_loop = 'loop-on-particles'
        pinCI.position_coordinates = ['x', 'y', 'z'] # determines the particle-storage dimensions
        pinCI.force_components = ['x', 'y','z']
        pinCI.force_precision = numpy.float64

        # Give the properties of the particle species.  The charges
        # and masses are normally those of the physical particles, and
        # not the computational macroparticles.  Macroparticle weights
        # are specified or computed in a separate file (see
        # user_particles_module below) giving the distribution
        # functions, and can vary from particle to particle.


        pinCI.particle_species = (('plasmaelectrons',
                             {'initial_distribution_type' : 'listed',
                              'charge' : -1.0*MyPlasmaUnits_C.elem_charge,
                              'mass' : 1.0*MyPlasmaUnits_C.electron_mass,
                              'dynamics' : 'implicit',
#                              'number_per_cell' : 12,
                              }
                             ),
                            ('Hplus', 
                             {'initial_distribution_type' : 'listed',
                              'charge' : 1.0*MyPlasmaUnits_C.elem_charge,
                              'mass' : 1.0*MyPlasmaUnits_C.AMU,
                              'dynamics' : 'implicit',
#                              'number_per_cell' : 6,
                              }
                             ),
                            ('He', 
                             {'initial_distribution_type' : None,
                              'charge' : 0.0,
                              'mass' : 4.0*MyPlasmaUnits_C.AMU,
                              'dynamics' : 'implicit',
#                              'number_per_cell' : 1,
                              }
                             ),
                            )

        # Give the name of the .py file containing the user-provided particle distributions
        pinCI.user_particles_module = "UserParticles_H_He_e"
#        print "setUp: UserParticle module is", pinCI.user_particles_module
        UPrt_M = im_M.import_module(pinCI.user_particles_module)
        pinCI.user_particles_class = UPrt_C = UPrt_M.ParticleDistributions_C

        self.pinCI = pinCI

#        self.particles = Particles_C(self.particle_species, self.pvars, precision, self.segment_length, UPrt_C, echoFlag=False)
        self.particles = Particle_C(pinCI, printFlag=False)

        return

    def test_species_names(self):
        """ Check that the species names are those specified by the user.
        """

        fncname = sys._getframe().f_code.co_name
        print '\ntest:', fncname, '('+__file__+')'

        particle_species = self.pinCI.particle_species
        # Check the names of the species
        for i in range(len(particle_species)):
            expected_name = particle_species[i][0]
            self.assertEqual(self.particles.species_names[i], expected_name, msg = "Species name is not correct")

    def test_listed_particles(self):
        """ 1. Check that the stored particles have the values listed by the user.
            2. Check the methods that count the particles.
        """

        fncname = sys._getframe().f_code.co_name
        print '\ntest:', fncname, '('+__file__+')'

        user_particles_class = self.pinCI.user_particles_class

        ninput_total = 0
        for sp in self.particles.species_names:
            init_dist_type = self.particles.initial_distribution_type[sp]
            if init_dist_type == 'listed':
                # Put user-listed particles into the storage array
                self.particles.create_from_list(sp, False)
                # Get the original input data by calling the
                # user-provided function in 'user_particles_file'.
                # The function name is the species name, and the
                # argument for the particles used in this particular
                # test is 'listed'
                ninput, particles = getattr(user_particles_class, sp)(init_dist_type)
                ninput_total += ninput
                # Check that the particles has the user-input values
                for ip in range(ninput):
                    putparticle = particles[ip]
                    # Get the stored particle data
#                    getparticle = self.particles.pseg_arr[sp][ip]
                    getparticle = self.particles.pseg_arr[sp].get(ip)
#                    print 'putparticle = ', putparticle
#                    print 'getparticle = ', getparticle
                    for ix in range(len(getparticle)):
                        self.assertAlmostEqual(getparticle[ix], putparticle[ix], msg="Particle is not correct")

                # Check the number of stored particles in each species
                nstored = self.particles.get_species_particle_count(sp, printFlag = False)
                self.assertEqual(nstored, ninput, msg = "Number of stored particles is not correct")

        # check total number of stored particles
        nstored = self.particles.get_total_particle_count(printFlag = False)
        self.assertEqual(nstored, ninput_total, msg = "Number of stored particles is not correct")

    def test_something(self):
        """ Check the number of particles in each species.
        """
        fncname = sys._getframe().f_code.co_name
        print '\ntest:', fncname, '('+__file__+')'

if __name__ == '__main__':
    unittest.main()
