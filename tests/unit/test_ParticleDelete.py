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

class TestParticleDeletion(unittest.TestCase):
    """Test deletion of particles from the storage arrays"""
    
    def setUp(self):

        # initializations for each test go here...

        '''
        self.segment_length = 100
        self.delete_flag = 0b01 # the lowest bit is 1
        self.trajectory_flag = 0b10 # the second lowest bit is 1

        precision = numpy.float64
        position_coordinates = ('x', 'y', 'z')
        momentum_coordinates = ('ux', 'uy', 'uz')
        phase_coordinates = position_coordinates + momentum_coordinates
        pvars = [coord for coord in phase_coordinates]

        pvars.append('weight')
        self.pvars = pvars
#        pvartypes = [precision for coord in phase_coordinates]
        pvartypes = [precision for var in pvars]

        pvars.append('bitflags')
        pvartypes.append(numpy.int32)

        pvars.append('cell_index')
        pvartypes.append(numpy.int32) # The size determines how many local cells you can have.

        self.particle_dtype = {'names' : pvars, 'formats': pvartypes}

        metadata = {'names' : pvars, 'formats': pvartypes}
        
        # Create the particle storage
        self.seg_array_obj = SegmentedArray_C(self.segment_length, metadata)

        # put a particle into the array

        '''


        # Create an instance of the DTparticleInput class
        pinCI = DTparticleInput_C()
        # Initialize particles
        pinCI.precision = numpy.float64
        pinCI.particle_integration_loop = 'loop-on-particles'

        pinCI.position_coordinates = ['x', 'y', 'z'] # determines the particle-storage dimensions

        pinCI.force_components = ['x', 'y', 'z']
        pinCI.force_precision = numpy.float64

        # Specify the particle species for this calculation

        # 1. electrons
        pinCI.particle_species = (('one_electron',
                             {'initial_distribution_type' : 'listed',
                              'charge' : -1.0*MyPlasmaUnits_C.elem_charge,
                              'mass' : 1.0*MyPlasmaUnits_C.electron_mass,
                              'dynamics' : 'explicit',
#                              'number_per_cell' : 12,
                              }
                             ),
        # 2. Hplus (proton)
                            ('Hplus', 
                             {'initial_distribution_type' : None,
                              'charge' : 1.0*MyPlasmaUnits_C.elem_charge,
                              'mass' : 1.0*MyPlasmaUnits_C.proton_mass,
                              'dynamics' : 'explicit',
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
        pinCI.user_particles_class = UPrt_M.UserParticleDistributions_C

        self.pinCI = pinCI

        # Make the storage array
        self.particleCI = Particle_C(pinCI, printFlag=False)

        # Initialize the particle properties
        for sp in self.particleCI.species_names:
            if self.particleCI.initial_distribution_type[sp] == 'listed':
                # Put user-listed particles into the storage array
                self.particleCI.create_from_list(sp, False)
        return
#    def setUp(self): #ENDDEF

    def test_1_delete_particles(self):
        """Delete some particles in the 'out' array"""

        fncname = sys._getframe().f_code.co_name
        print '\ntest: ', fncname, '('+__file__+')'

        ctrlCI = DT_M.DTcontrol_C()

        ctrlCI.dt = 1.0e-5
        ctrlCI.nsteps = 1

        # Set constant fields
        E0 = (1.0e-4, 2.0e-4, 3.0e-4)
#        B0 = (0.0, 0.0, 0.0)
        ctrlCI.E0 = Vec_C(self.particleCI.dimension, E0)
#        ctrlCI.B0 = Vec_C(self.particleCI.dimension, B0)

        # Create more particles to test deletion. Copy the particle
        # that's already stored to make 3 full and one partially
        # filled segments.
        # We end up with 306 particles, requiring 4 segments
        num_particles = 5+3*self.particleCI.SEGMENT_LENGTH
        dx = 0.2
        for sp in self.particleCI.species_names:
            if self.particleCI.get_species_particle_count(sp) == 0: continue # Skip if there are no particles in this species
            getparticle = self.particleCI.pseg_arr[sp].get(0) # why get(0) instead of [0]?
            putparticle = getparticle
            x = putparticle[0]
            for i in range(num_particles):
#                print 'putparticle = ', putparticle
                putparticle[0] = x
    #            putparticle = (x, y, z, px, py, pz, weight, bitflags, cell_index,)
                self.particleCI.pseg_arr[sp].put(putparticle)
                x += dx

        # Query the particle arrays after initialization.
        for sp in self.particleCI.species_names:
            (nseg_in, nseg_out) = self.particleCI.pseg_arr[sp].get_number_of_segments()
            print sp, "species has %d segments in the 'in' array and %d segments in the 'out' array" % (nseg_in, nseg_out)
#            npart_out = self.particleCI.pseg_arr[sp].get_number_of_items()
            npart_out = self.particleCI.get_species_particle_count(sp)
            print sp, "species has %d particles in the 'out' array" % npart_out

        # Delete some particles

        for sp in self.particleCI.species_names:
            np = self.particleCI.get_species_particle_count(sp)
            if np == 0: continue # Skip if there are no particles in this species
            # Delete 7 particles, which fit into 3 segments
            delparts = (0, np-1, np/51, np/31, np/11, np/5, np/3)
            # Delete 6 particles, which fit into 3 segments
#            delparts = (0, np-1, np/51, np/31, np/11, np/5)
            # Delete 5 particles, which needs a 4th segment
#            delparts = (0, np-1, np/51, np/31, np/11)
            for ip in delparts:
                getparticle = self.particleCI.pseg_arr[sp].get(ip) # why get(0) instead of [0]? more explicit?
#                print "getparticle['bitflags']", getparticle['bitflags']
                # Set the delete flag for these particles.
                # Since getparticle is a reference, this modifies the
                # stored particle.
                getparticle['bitflags'] = getparticle['bitflags'] | Particle_C.DELETE_FLAG
#                print "getparticle['bitflags']", getparticle['bitflags']

        # Move the particles one timestep. This serves to copy the
        # particles from 'in' to the 'out' array.  Don't care about
        # the actual motion

        ncoords = self.particleCI.dimension # number of particle coordinates to check
#        isp = 0
        print "Moving", self.particleCI.get_total_particle_count(), "particles for", ctrlCI.nsteps, "timesteps"
        for sp in self.particleCI.species_names:
            if self.particleCI.get_species_particle_count(sp) == 0: continue # Skip if there are no particles in this species
            self.particleCI.move_particles_in_uniform_fields(sp, ctrlCI)
            # Check that the first particles in the array reach the right speed
#            getparticle = self.particleCI.pseg_arr[sp][0]
            getparticle = self.particleCI.pseg_arr[sp].get(0)
#            print 'calculated = ', getparticle
#            print 'expected = ', p_expected[isp]
            for ic in range(ncoords):
#            for ix in range(len(getparticle)):
#                self.assertAlmostEqual(getparticle[ic], p_expected[isp][ic], msg="Particle is not in correct position")
                pass
#            isp += 1


        # Now check on the arrays again
        print "\nAfter a particle move step:\n"
        for sp in self.particleCI.species_names:
            (nseg_in, nseg_out) = self.particleCI.pseg_arr[sp].get_number_of_segments()
            print sp, "species has %d segments in the 'in' array and %d segments in the 'out' array" % (nseg_in, nseg_out)
#            npart_out = self.particleCI.pseg_arr[sp].get_number_of_items()
            npart_out = self.particleCI.get_species_particle_count(sp)
            print sp, "species has %d particles in the 'out' array" % npart_out

        return
#    def test_1_delete_particles(self): ENDDEF

#class TestParticleDeletion(unittest.TestCase): ENDCLASS

if __name__ == '__main__':
    unittest.main()
