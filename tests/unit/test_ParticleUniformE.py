#!/usr/bin/env python

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import numpy
import importlib as im_M
import unittest

import DT_Module as DT_M

from DT_Module import DTparticleInput_C
from UserUnits_Module import MyPlasmaUnits_C
from Particle_Module import Particle_C

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

class TestParticleUniformE(unittest.TestCase):
    """Test classes in Particle_Module"""
    
    def setUp(self):

        # initializations for each test go here...

        # Create an instance of the DTparticleInput class
        pinCI = DTparticleInput_C()
        # Initialize particles
        pinCI.precision = numpy.float64
        pinCI.particle_integration_loop = 'loop-on-particles'
        pinCI.position_coordinates = ['x', 'y', 'z'] # determines the particle-storage dimensions
        pinCI.force_components = ['x', 'y',]
        pinCI.force_precision = numpy.float64

        # Specify the particle species for this calculation
        # 1. electrons
        pinCI.particle_species = (('plasmaelectrons',
                             {'initial_distribution_type' : 'listed',
                              'charge' : -1.0*MyPlasmaUnits_C.elem_charge,
                              'mass' : 1.0*MyPlasmaUnits_C.electron_mass,
                              'dynamics' : 'explicit',
#                              'number_per_cell' : 12,
                              }
                             ),
        # 2. Hplus (proton)
                            ('Hplus', 
                             {'initial_distribution_type' : 'listed',
                              'charge' : 1.0*MyPlasmaUnits_C.elem_charge,
                              'mass' : 1.0*MyPlasmaUnits_C.proton_mass,
                              'dynamics' : 'explicit',
#                              'mass' : 1.0*MyPlasmaUnits_C.AMU,
#                              'number_per_cell' : 6,
                              }
                             ),
        # 3. Neutral: test when there are no particles
                            ('He', 
                             {'initial_distribution_type' : None,
                              'charge' : 0.0,
                              'mass' : 4.0*MyPlasmaUnits_C.AMU,
                              'dynamics' : 'explicit',
#                              'number_per_cell' : 1,
                              }
                             ),
                            )

        # Provide the initial conditions for the above species
        pinCI.user_particles_module = "UserParticles_H_He_e"
        UPrt_M = im_M.import_module(pinCI.user_particles_module)
        pinCI.user_particles_class = UPrt_M.ParticleDistributions_C

        self.pinCI = pinCI

        # Make the storage array
        self.particleCI = Particle_C(pinCI, printFlag=False)

        # Store the particles
        for sp in self.particleCI.species_names:
            if self.particleCI.initial_distribution_type[sp] == 'listed':
                # Put user-listed particles into the storage array
                self.particleCI.create_from_list(sp, False)
        return

    def test_1_electric_field_push_1step(self):
        """ Check that the electric field push is correct.  Push
            sample particles for 1 step.
        """
        fncname = sys._getframe().f_code.co_name
        print '\ntest: ', fncname

        ctrlCI = DT_M.DTcontrol_C()

        ctrlCI.dt = 1.0e-5
        ctrlCI.nsteps = 1

        # Set constant fields
        E0 = (1.0e-4, 2.0e-4, 3.0e-4)
#        B0 = (0.0, 0.0, 0.0)
        ctrlCI.E0 = Vec_C(self.particleCI.dimension, E0)
#        ctrlCI.B0 = Vec_C(self.particleCI.dimension, B0)

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
        ncoords = self.particleCI.dimension # number of particle coordinates to check
        isp = 0
        for sp in self.particleCI.species_names:
            if self.particleCI.get_species_particle_count(sp) == 0: continue
            self.particleCI.move_particles_in_uniform_fields(sp, ctrlCI)
            # Check that the first particles in the array reach the right speed
#            getparticle = self.particleCI.pseg_arr[sp][0]
            getparticle = self.particleCI.pseg_arr[sp].get(0)
#            print 'calculated = ', getparticle
#            print 'expected = ', p_expected[isp]
            for ic in range(ncoords):
#            for ix in range(len(getparticle)):
                self.assertAlmostEqual(getparticle[ic], p_expected[isp][ic], msg="Particle is not in correct position")
            isp += 1

        return

    def test_2_electric_field_push_10steps(self):
        """ Check that the electric field push is correct.  Push
            sample particles for 10 steps.
        """
        fncname = sys._getframe().f_code.co_name
        print '\ntest: ', fncname

        ctrlCI = DT_M.DTcontrol_C()

        ctrlCI.dt = 1.0e-5
        ctrlCI.nsteps = 10

        # Set constant fields
        E0 = (1.0e-4, 2.0e-4, 3.0e-4)
#        B0 = (0.0, 0.0, 0.0)
        ctrlCI.E0 = Vec_C(self.particleCI.dimension, E0)
#        ctrlCI.B0 = Vec_C(self.particleCI.dimension, B0)

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

        # Integrate for nsteps
        for istep in xrange(ctrlCI.nsteps):
#            print 'test_2_electric_field_push_10steps: istep =', istep
            # Loop on the species
            for sp in self.particleCI.species_names:
#                if sp != 'plasmaelectrons': continue
#                print 'test_2_electric_field_push_10steps: sp =', sp
                if self.particleCI.get_species_particle_count(sp) == 0: continue
                self.particleCI.move_particles_in_uniform_fields(sp, ctrlCI)

        # check the final position of the first particle of each species
        ncoords = self.particleCI.dimension # number of particle coordinates to check
        isp = 0
        for sp in self.particleCI.species_names:
            if self.particleCI.get_species_particle_count(sp) == 0: continue
            # Check that the first particles in the array reach the right speed
            getparticle = self.particleCI.pseg_arr[sp].get(0)
#            print 'calculated = ', getparticle
#            print 'expected = ', p_expected[isp]
            for ic in range(ncoords):
#            for ix in range(len(getparticle)):
                self.assertAlmostEqual(getparticle[ic], p_expected[isp][ic], msg="Particle is not in correct position")
            isp += 1

        return

    def test_3_magnetic_field_push(self):            
        """ 1. Check that the magnetic field push is correct.
            2. 
        """
        fncname = sys._getframe().f_code.co_name
        print '\ntest: ', fncname
        pass

        return

    def test_4_something(self):
        """ Check the number of particles in each species.
        """
        fncname = sys._getframe().f_code.co_name
        print '\ntest: ', fncname
        pass

        return
# class TestParticleUniformE(unittest.TestCase): ENDCLASS


if __name__ == '__main__':
    unittest.main()
