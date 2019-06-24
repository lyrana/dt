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
from RecordedData_Module import *

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
class TestParticleDeletion(unittest.TestCase):
    """Test deletion of particles from the storage arrays"""
    
    def setUp(self):

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

        # initializations for each test go here...

        self.ctrl = DT_m.DTcontrol_C()

        self.ctrl.time = 0.0
        self.ctrl.dt = 1.0e-5
        self.ctrl.n_timesteps = 1

        ### Create an instance of the DTparticleInput class

        pin = ParticleInput_C()

        # Set particle parameters
        pin.precision = numpy.float64
        pin.particle_integration_loop = 'loop-on-particles'

        pin.coordinate_system = 'cartesian_x_y_z'
# TODO: Set these based on the coordinate_system:
        pin.position_coordinates = ['x', 'y', 'z'] # determines the particle-storage dimensions

        pin.force_components = ['x', 'y', 'z']
        pin.force_precision = numpy.float64

        ### Specify the particle species for this calculation

        # 1. electrons
#         pin.particle_species = (('one_electron',
#                              {'initial_distribution_type' : 'listed',
#                               'charge' : -1.0*MyPlasmaUnits_C.elem_charge,
#                               'mass' : 1.0*MyPlasmaUnits_C.electron_mass,
#                               'dynamics' : 'explicit',
# #                              'number_per_cell' : 12,
#                               }
#                              ),
#         # 2. Hplus (proton)
#                             ('Hplus', 
#                              {'initial_distribution_type' : None,
#                               'charge' : 1.0*MyPlasmaUnits_C.elem_charge,
#                               'mass' : 1.0*MyPlasmaUnits_C.proton_mass,
#                               'dynamics' : 'explicit',
# #                              'number_per_cell' : 6,
#                               }
#                              ),
#         # 3. Neutral: test when there are no particles
#                             ('He', 
#                              {'initial_distribution_type' : None,
#                               'charge' : 0.0,
#                               'mass' : 4.0*MyPlasmaUnits_C.AMU,
#                               'dynamics' : 'explicit',
# #                              'number_per_cell' : 1,
#                               }
#                              ),

#                             )

        speciesName = 'one_electron'
        charge = -1.0*MyPlasmaUnits_C.elem_charge
        mass = 1.0*MyPlasmaUnits_C.electron_mass
        dynamics = 'explicit'
        oneElectron_S = ParticleSpecies_C(speciesName, charge, mass, dynamics)

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
        pin.particle_species = (oneElectron_S, HPlus_S, He_S,
                                 )
        # Make the particle object from pin
        self.particle_P = Particle_C(pin, print_flag=False)

        # Give the name of the .py file containing additional particle data (lists of
        # particles, boundary conditions, source regions, etc.)
        userParticlesModuleName = "UserParticles_H_He_e"

        # Import this module
        userParticlesModule = im_m.import_module(userParticlesModuleName)

        # self.particle_P.user_particles_module_name = userParticlesModuleName
        self.particle_P.user_particles_class = userParticlesClass = userParticlesModule.UserParticleDistributions_C

        ### The one_electron species is present at t=0. The other two species are not.

        # Name the initialized species (it should be in species_names above)
        speciesName = 'one_electron'
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
        oneElectronParams = {'species_name': speciesName,
                              'initial_distribution_type': initialDistributionType,
                              }

        # The dictionary keys are mnemonics for the initialized particles
        initialParticlesDict = {'initial_one_electron': (oneElectronParams,),
                                }

        self.particle_P.initial_particles_dict = initialParticlesDict


        ### Create a particle trajectory object

        # Use an input object to collect initialization data for the trajectory object
        self.trajin = TrajectoryInput_C()

        self.trajin.maxpoints = None # Set to None to get every point
        self.trajin.extra_points = 1  # Set to 1 to make sure one boundary-crossing can be
                                      # accommodated. Set to a larger value if there are
                                      # multiple boundary reflections.

        # Specify which particle variables to save.  This has the
        # form of a numpy dtype specification.
        
        name_list_base = ['step', 't'] # These are always recorded
        format_list_base = [int, numpy.float32] # Start off the format list with types for
                                                # 'step' and 't'

        explicit_attributes = ['x', 'ux', 'y', 'uy', 'Ex', 'Ey']
        format_list = format_list_base + [numpy.float32]*len(explicit_attributes)
        self.trajin.explicit_dict = {'names': name_list_base+explicit_attributes, 'formats': format_list}
        
        implicit_attributes = ['x', 'ux', 'phi']
        format_list = format_list_base + [numpy.float32]*len(implicit_attributes)
        self.trajin.implicit_dict = {'names': name_list_base+implicit_attributes, 'formats': format_list}

        neutral_attributes = ['x', 'ux', 'y', 'uy']
        format_list = format_list_base + [numpy.float32]*len(neutral_attributes)
        self.trajin.neutral_dict = {'names': name_list_base+neutral_attributes, 'formats': format_list}

        ## Create the trajectory object and attach it to the particle object.
        # No trajectory storage is created until particles
        # with TRAJECTORY_FLAG on are encountered.
        p_P = self.particle_P # abbreviation
        traj_T = Trajectory_C(self.trajin, self.ctrl, p_P.explicit_species, p_P.implicit_species, p_P.neutral_species)
        self.particle_P.traj_T = traj_T

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
#    def setUp(self): #ENDDEF

    def test_1_delete_particles(self):
        """Delete some particles in the 'out' array"""

        fncName = sys._getframe().f_code.co_name
        print('\ntest: ', fncName, '('+__file__+')')

        # Abbreviations

        p_P = self.particle_P
        traj_T = p_P.traj_T

        # Set constant fields
        E0 = (1.0e-4, 2.0e-4, 3.0e-4)
#        B0 = (0.0, 0.0, 0.0)
        self.ctrl.E0 = Vec_C(p_P.particle_dimension, E0)
#        ctrl.B0 = Vec_C(p_P.particle_dimension, B0)

        print("SEGMENT_LENGTH is %d" % p_P.SEGMENT_LENGTH)

        # Create more particles to test deletion. Copy the particle
        # that's already stored to make 3 full and one partially
        # filled segments.
        # We end up with 306 particles, requiring 4 segments
        num_particles = 5+3*p_P.SEGMENT_LENGTH
        dx = 0.2
        for sp in p_P.species_names:
            if p_P.get_species_particle_count(sp) == 0: continue # Skip if there are no particles in this species
# If we use get_item(), we get a dict, which doesn't have the right type for the following.
#            getparticle = p_P.pseg_arr[sp].get_item(0)

            getparticle = p_P.pseg_arr[sp].get_as_list(0)
            
#            (pseg, offset) = p_P.pseg_arr[sp].get_segment_and_offset(0)
#            getparticle = pseg[offset]
            
            print ("getparticle = ", getparticle)
            putparticle = getparticle
            x = putparticle[0]
            print("Add %d more particles to %s species" % (num_particles, sp))
            for i in range(num_particles):
#                print 'putparticle = ', putparticle
                # Here's the layout of data in the putparticle tuple:
                #   (x, y, z, x0, y0, z0, px, py, pz, weight, bitflags, cell_index,)
                putparticle[0] = x
                # Turn on trajectory flag for a particle near the end:
                if i == num_particles-10:
                    putparticle[10] = 0b00 | p_P.TRAJECTORY_FLAG # turn on trajectory flag
                else:
                    putparticle[10] = 0b00 # initialize all bits to 0

                # Store the particle
#                p, pindex = p_P.pseg_arr[sp].push_back(putparticle)
                seg_index, full_index = p_P.pseg_arr[sp].push_back(putparticle)
#                p = p_P.pseg_arr[sp][full_index] # Retrieve the particle from the SAP.

                # Retrieve the particle structure using the returned Numpy array it's in.
                (parray, offset) = p_P.pseg_arr[sp].get_segment_and_offset(full_index)
                p = parray[offset] # Retrieve the particle from the SAP.
                
                ## Add the new particles to the trajectory object.
                if p['bitflags'] & p_P.TRAJECTORY_FLAG != 0:
                    if traj_T is not None:
                        print('full_index for trajectory = ', full_index)
                        traj_T.particle_index_list[sp].append(full_index)
                        dynamicsType = 'explicit'
                        traj_T.create_trajectory(sp, full_index, dynamicsType)
                    else:
    # Instead of printing this message, a traj_T object could be created here.
                        print(fncName, "*** DT Warning: A trajectory flag is on, but no trajectory object has been created yet. ***")
                x += dx

        # Query the particle arrays after initialization.
        for sp in p_P.species_names:
            (nseg_in, nseg_out) = p_P.pseg_arr[sp].get_number_of_segments()
            print(sp, "species has %d segments in the 'in' array and %d segments in the 'out' array" % (nseg_in, nseg_out))
#            npart_out = p_P.pseg_arr[sp].get_number_of_items()
            npart_out = p_P.get_species_particle_count(sp)
            print(sp, "species has %d particles in the 'out' array" % npart_out)

        # Delete some of the particles in each species.

        for sp in p_P.species_names:
            np = p_P.get_species_particle_count(sp)
            if np == 0: continue # Skip if there are no particles in this species
            # Delete 7 particles, which fit into 3 segments
#26aug18: In Python3, the numbers are upgraded to float before dividing
#            delparts = (0, np-1, np/51, np/31, np/11, np/5, np/3) # This is (0, 0, ...)
                                                                   # for np=1
#            print("delparts at a =", delparts)
            delparts = (0, np-1, int(np/51), int(np/31), int(np/11), int(np/5), int(np/3))
#            print("delparts at b =", delparts)
            
            # Delete 6 particles, which fit into 3 segments
#            delparts = (0, np-1, np/51, np/31, np/11, np/5)
            # Delete 5 particles, which needs a 4th segment
#            delparts = (0, np-1, np/51, np/31, np/11)
            for ip in delparts:
                # We need to change the stored data, so use get_segment_and_offset()
                (parray, offset) = p_P.pseg_arr[sp].get_segment_and_offset(full_index)
                getparticle = parray[offset]
# This systax only works for the Python SAP:                
#                getparticle = p_P.pseg_arr[sp].get_item(ip)

#                print ("getparticle['bitflags']", getparticle['bitflags'])
                # Set the delete flag for these particles.
                # Since getparticle returns a REFERENCE, this modifies the
                # stored particle.
                getparticle['bitflags'] = getparticle['bitflags'] | Particle_C.DELETE_FLAG
#                print ("getparticle['bitflags']", getparticle['bitflags'])

        # Move the particles one timestep. This serves to copy the particles from 'in' to
        # the 'out' array.  The deleted particles will get dropped from the 'out' array.
        # Don't care about the actual motion.

        ncoords = p_P.particle_dimension # number of particle coordinates to check
#        isp = 0
        print("Moving", p_P.get_total_particle_count(), "particles for", self.ctrl.n_timesteps, "timesteps")
        for sp in p_P.species_names:
            if p_P.get_species_particle_count(sp) == 0: continue # Skip if there are no particles in this species
            p_P.move_particles_in_uniform_fields(sp, self.ctrl)

        # Now check on the arrays again
        print("\nAfter a particle move step:\n")
        for sp in p_P.species_names:
            (nseg_in, nseg_out) = p_P.pseg_arr[sp].get_number_of_segments()
            print(sp, "species has %d segments in the 'in' array and %d segments in the 'out' array" % (nseg_in, nseg_out))
            npart_out = p_P.get_species_particle_count(sp)
            print(sp, "species has %d particles in the 'out' array" % npart_out)

# Need an assert test here.

        return
#    def test_1_delete_particles(self): ENDDEF

#class TestParticleDeletion(unittest.TestCase): ENDCLASS

if __name__ == '__main__':
    unittest.main()
