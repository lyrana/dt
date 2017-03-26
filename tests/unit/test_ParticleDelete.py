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
from DT_Module import DTtrajectoryInput_C

from UserUnits_Module import MyPlasmaUnits_C
from Particle_Module import Particle_C
from Trajectory_Module import Trajectory_C

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


        self.ctrlCI = DT_M.DTcontrol_C()

        self.ctrlCI.dt = 1.0e-5
        self.ctrlCI.nsteps = 1

        ### Create an instance of the DTparticleInput class

        pinCI = DTparticleInput_C()

        # Set particle parameters
        pinCI.precision = numpy.float64
        pinCI.particle_integration_loop = 'loop-on-particles'

        pinCI.position_coordinates = ['x', 'y', 'z'] # determines the particle-storage dimensions

        pinCI.force_components = ['x', 'y', 'z']
        pinCI.force_precision = numpy.float64

        ### Specify the particle species for this calculation

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

        ### Create a particle trajectory object

        # Use an input object to collect initialization data for the trajectory object
        self.trajinCI = DTtrajectoryInput_C()

        self.trajinCI.maxpoints = None # Set to None to get every point

        # Specify which particle variables to save.  This has the
        # form of a numpy dtype specification.
        self.trajinCI.explicitDict = {'names': ['x', 'ux', 'y', 'uy', 'Ex', 'Ey'], 'formats': [numpy.float32]*6}
        self.trajinCI.implicitDict = {'names': ['x', 'ux', 'phi'], 'formats': [numpy.float32]*3}
        self.trajinCI.neutralDict = {'names': ['x', 'ux', 'y', 'uy'], 'formats': [numpy.float32]*4}

        ## Create the trajectory object and attach it to the particle object.
        # No trajectory storage is created until particles
        # with TRAJECTORY_FLAG on are encountered.
        pCI = self.particleCI # abbreviation
        trajCI = Trajectory_C(self.trajinCI, self.ctrlCI, pCI.explicit_species, pCI.implicit_species, pCI.neutral_species)
        self.particleCI.trajCI = trajCI

        # Initialize the particle properties
        for sp in self.particleCI.species_names:
            if self.particleCI.initial_distribution_type[sp] == 'listed':
                # Put user-listed particles into the storage array
                self.particleCI.create_from_list(sp, False)
        return
#    def setUp(self): #ENDDEF

    def test_1_delete_particles(self):
        """Delete some particles in the 'out' array"""

        fncName = sys._getframe().f_code.co_name
        print '\ntest: ', fncName, '('+__file__+')'

        # Abbreviations

        pCI = self.particleCI
        trajCI = pCI.trajCI

        # Set constant fields
        E0 = (1.0e-4, 2.0e-4, 3.0e-4)
#        B0 = (0.0, 0.0, 0.0)
        self.ctrlCI.E0 = Vec_C(pCI.particle_dimension, E0)
#        ctrlCI.B0 = Vec_C(pCI.particle_dimension, B0)

        print "SEGMENT_LENGTH is %d" % pCI.SEGMENT_LENGTH

        # Create more particles to test deletion. Copy the particle
        # that's already stored to make 3 full and one partially
        # filled segments.
        # We end up with 306 particles, requiring 4 segments
        num_particles = 5+3*pCI.SEGMENT_LENGTH
        dx = 0.2
        for sp in pCI.species_names:
            if pCI.get_species_particle_count(sp) == 0: continue # Skip if there are no particles in this species
            getparticle = pCI.pseg_arr[sp].get(0) # why get(0) instead of [0]?
            putparticle = getparticle
            x = putparticle[0]
            print "Add %d more particles to %s species" % (num_particles, sp)
            for i in range(num_particles):
#                print 'putparticle = ', putparticle
                # Here's the layout of data in the putparticle tuple:
                #   (x, y, z, x0, y0, z0, px, py, pz, weight, bitflags, cell_index,)
                putparticle[0] = x
                # Turn on trajectory flag for a particle near the end:
                if i == num_particles-10:
                    putparticle[10] = 0b00 | pCI.TRAJECTORY_FLAG # turn on trajectory flag
                else:
                    putparticle[10] = 0b00 # initialize all bits to 0

                # Store the particle
                p, pindex = pCI.pseg_arr[sp].put(putparticle)

                ## Add the new particles to the trajectory object.
                if p['bitflags'] & pCI.TRAJECTORY_FLAG != 0:
                    if trajCI is not None:
                        print 'pindex for trajectory = ', pindex
                        trajCI.ParticleIdList[sp].append(pindex)
                        dynamics_type = 'explicit'
                        trajCI.create_trajectory(sp, dynamics_type)
                    else:
    # Instead of printing this message, a trajCI object could be created here.
                        print fncName, "*** DT Warning: A trajectory flag is on, but no trajectory object has been created yet. ***"
                x += dx

        # Query the particle arrays after initialization.
        for sp in pCI.species_names:
            (nseg_in, nseg_out) = pCI.pseg_arr[sp].get_number_of_segments()
            print sp, "species has %d segments in the 'in' array and %d segments in the 'out' array" % (nseg_in, nseg_out)
#            npart_out = pCI.pseg_arr[sp].get_number_of_items()
            npart_out = pCI.get_species_particle_count(sp)
            print sp, "species has %d particles in the 'out' array" % npart_out

        # Delete some of the particles in each species.

        for sp in pCI.species_names:
            np = pCI.get_species_particle_count(sp)
            if np == 0: continue # Skip if there are no particles in this species
            # Delete 7 particles, which fit into 3 segments
            delparts = (0, np-1, np/51, np/31, np/11, np/5, np/3) # This is (0, 0, ...)
                                                                  # for np=1
            # Delete 6 particles, which fit into 3 segments
#            delparts = (0, np-1, np/51, np/31, np/11, np/5)
            # Delete 5 particles, which needs a 4th segment
#            delparts = (0, np-1, np/51, np/31, np/11)
            for ip in delparts:
                getparticle = pCI.pseg_arr[sp].get(ip) # why get(0) instead of [0]? more explicit?
#                print "getparticle['bitflags']", getparticle['bitflags']
                # Set the delete flag for these particles.
                # Since getparticle returns a REFERENCE, this modifies the
                # stored particle.
                getparticle['bitflags'] = getparticle['bitflags'] | Particle_C.DELETE_FLAG
#                print "getparticle['bitflags']", getparticle['bitflags']

        # Move the particles one timestep. This serves to copy the particles from 'in' to
        # the 'out' array.  The deleted particles will get dropped from the 'out' array.
        # Don't care about the actual motion.

        ncoords = pCI.particle_dimension # number of particle coordinates to check
#        isp = 0
        print "Moving", pCI.get_total_particle_count(), "particles for", self.ctrlCI.nsteps, "timesteps"
        for sp in pCI.species_names:
            if pCI.get_species_particle_count(sp) == 0: continue # Skip if there are no particles in this species
            pCI.move_particles_in_uniform_fields(sp, self.ctrlCI)

        # Now check on the arrays again
        print "\nAfter a particle move step:\n"
        for sp in pCI.species_names:
            (nseg_in, nseg_out) = pCI.pseg_arr[sp].get_number_of_segments()
            print sp, "species has %d segments in the 'in' array and %d segments in the 'out' array" % (nseg_in, nseg_out)
            npart_out = pCI.get_species_particle_count(sp)
            print sp, "species has %d particles in the 'out' array" % npart_out

# Need an assert test here.

        return
#    def test_1_delete_particles(self): ENDDEF

#class TestParticleDeletion(unittest.TestCase): ENDCLASS

if __name__ == '__main__':
    unittest.main()
