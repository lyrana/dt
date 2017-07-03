#!/usr/bin/env python

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import os
import importlib as im_M
import unittest

import numpy as np_M

import dolfin as df_M

from DT_Module import DTcontrol_C

from Dolfin_Module import *
from Particle_Module import *
from UserMesh_FE_XYZ_Module import *

from UserUnits_Module import MyPlasmaUnits_C

#STARTCLASS
class TestParticleInitialization(unittest.TestCase):
    """Test classes in Particle_Module"""
    
    def setUp(self):

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

        # Initializations performed before each test go here...

        # Create an instance of the DTparticleInput class
        pinCI = self.pinCI = ParticleInput_C()
        # Set up particle variables
        pinCI.precision = np_M.float64

        pinCI.particle_integration_loop = 'loop-on-particles'
        pinCI.position_coordinates = ['x', 'y', 'z'] # determines the particle-storage dimensions
        pinCI.force_components = ['x', 'y', 'z']
        pinCI.force_precision = np_M.float64


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
        pCI = self.particleCI = Particle_C(pinCI, print_flag=False)

        # Give the name of the .py file containing special particle data (lists
        # of particles, boundary conditions, particle-initialization regions,
        # particle source regions, etc.)
        userParticleModule = "UserParticles_H_He_e"

        # Import this module
        UPrt_M = im_M.import_module(userParticleModule)

        pCI.user_particle_module = userParticleModule
        pCI.user_particle_class = userParticleClass = UPrt_M.UserParticleDistributions_C

        ### Provide input for particles present at t=0

        ## plasma electrons are present at t=0
        # Name the initialized species (it should be in species_names above)
        speciesName = 'plasma_electrons'
        # Check that this species has been defined above
        if speciesName not in pCI.species_names:
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
        if speciesName not in pCI.species_names:
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

        pCI.initial_particles_dict = initialParticlesDict

        return
#    def setUp(self):ENDDEF

#class TestParticleInitialization(unittest.TestCase):
    def test_species_names(self):
        """ Check that the species names are those specified by the user.
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest:', fncName, '('+__file__+')'

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

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest: ', fncName, '('+__file__+')'

        userParticleClass = self.particleCI.user_particle_class

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
                # Get the original listed-particle input data by calling the
                # user-provided function in 'user_particle_file'.
                # The function name is the species name, and the
                # argument for the particles used in this particular
                # test is 'listed'
                ninput, particles = getattr(userParticleClass, s)(initialDistributionType)
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
                nStored = self.particleCI.get_species_particle_count(s, print_flag = False)
                self.assertEqual(nStored, ninput, msg = "Number of stored particles is not correct")

        # check total number of stored particles
        nStored = self.particleCI.get_total_particle_count(print_flag = False)
        self.assertEqual(nStored, ninput_total, msg = "Number of stored particles is not correct")

        return
#    def test_listed_particles(self):ENDDEF

#class TestParticleInitialization(unittest.TestCase):
    def test_functional_particles(self):
        """ Test initialization of particles using a function over a region of
            the mesh.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest: ', fncName, '('+__file__+')'

        ## Control struct
        ctrlCI = DTcontrol_C()

        ## Mesh input
        uMI2DCI = UserMeshInput_C()

# Replace with TUPLES since we're at toplevel
        uMI2DCI.pmin = df_M.Point(-10.0, -10.0)
        uMI2DCI.pmax = df_M.Point(10.0, 10.0)
        uMI2DCI.cells_on_side = (2, 2)
        uMI2DCI.diagonal = 'left'

        ## Create a 2D particle mesh, since we're going to initialize particles
        ## on a region of a mesh.
        # UserMesh_FE_XYZ_Module can make the mesh from the above input.
        plotFlag = False
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name + ": mesh"
        pmesh2DCI = UserMesh_C(uMI2DCI, compute_dictionaries=True, compute_tree=True, plot_flag=plotFlag, plot_title=plotTitle)

        pCI = self.particleCI

        userParticleClass = pCI.user_particle_class

        ## 1. Initial hot electrons on bottom-left side of the mesh

        # a. Provide the species name and physical parameters
        #    for this particle-initialization
        speciesName = 'plasma_electrons'

        # Check that this species has been defined
        if speciesName in pCI.species_names:
            charge = pCI.charge[speciesName]
            mass = pCI.mass[speciesName]
        else:
            print "The species", speciesName, "needs to be defined."
            sys.exit()

        initialDistributionType = 'function_over_region'
        # Specify the initial physical number-density for this species
        numberDensity = 1.0e12
        # Check for positivity
        if numberDensity <= 0:
            print "Check number density for species", speciesName, "is zero or negative. Should be positive"
            sys.exit()

        # Compute a value for thermalSpeed from the temperature in eV
        temperature_eV = 2.0
        temp_joule = temperature_eV*MyPlasmaUnits_C.elem_charge
        thermalSpeed = np_M.sqrt(2.0*temp_joule/mass)

        # Set a drift velocity
        # Set a drift velocity
        vdrift_x = 1.0
        vdrift_y = 2.0
        vdrift_z = 3.0

        driftVelocity = (vdrift_x, vdrift_y, vdrift_z)
#        driftVelocity[:] = (1.0,)

        # Set desired particles-per-cell
        numberPerCell = 1000

        # Collect the parameters for this initialization in a dictionary
        hotElectronParams = {'species_name': speciesName,
                             'initial_distribution_type': initialDistributionType,
                             'number_density': numberDensity,
                             'thermal_speed': thermalSpeed,
                             'drift_velocity': driftVelocity,
                             'number_per_cell': numberPerCell}

        # b. Name the particle-creation function.
        maxwellianGenerator = pCI.create_maxwellian_particles

        # c. Initialization-region geometry
        xmin = -9.0; ymin = -9.0
        xmax = -4.0; ymax = -4.0

# Don't use df_M in the toplevel routine, if possible
#        pmin = df_M.Point(xmin, ymin)
#        pmax = df_M.Point(xmax, ymax)
        pmin = (xmin, ymin)
        pmax = (xmax, ymax)
        hotElectronsRegionCI = RectangularRegion_C(pmesh2DCI, pmin, pmax)

        ## Put all the particle-initialization data into a dictionary

        # The dictionary keys are mnemonics for particle-initialization names.
        # The dictionary values are lists of tuples giving the 
        #   a. the physical initialization parameters.
        #   b. particle creation function and
        #   c. initialization region

        ## Note: the names here are particle-initialization names, NOT species names!
        initialParticlesDict = {'initial_hot_electrons': (hotElectronParams, maxwellianGenerator, hotElectronsRegionCI),
#                              'initial_background_electrons': (backgroundElectronParams, maxwellianGenerator, wholeMesh),
                              }

        # Replace the 'listed' initial particles made in Setup() with
        # the 'functional' initial particles created above

        pCI.initial_particles_dict = initialParticlesDict

        # Loop on the initialization methods in the dictionary

        # For 'listed' initialization, the dictionary has entries like
        # {ipName: (ipParams,)}
        # For 'function_over_region' initialization, the dictionary has entries
        # like
        # {ipName: (ipParams, ipFunc, ipRegionCI)}

        nCreatedTotal = 0
        for ipName, ipTuple in initialParticlesDict.iteritems():
            ipParams = ipTuple[0]
            s = ipParams['species_name']
            initialDistributionType = ipParams['initial_distribution_type']

            if initialDistributionType == 'function_over_region':
                print fncName, "Initializating", ipName, "particles"
                ipFunc = ipTuple[1]
                ipRegionCI = ipTuple[2]
                # Invoke the creation function
                time = 0.0
                ipFunc(time, ipRegionCI, ipParams)

                # Check the number of stored particles for each species

                nCreated = ipRegionCI.ncell * ipParams['number_per_cell']
                print fncName, "Number of particles created =", nCreated
                nCreatedTotal += nCreated

                nStored = pCI.get_species_particle_count(s, print_flag = False)
                self.assertEqual(nStored, nCreated, msg = "Number of stored particles is not correct")

        # check total number of stored particles
        nStored = pCI.get_total_particle_count(print_flag = False)
        print fncName, "Total number of particles created =", nCreatedTotal
        self.assertEqual(nStored, nCreatedTotal, msg = "Total number of stored particles is not correct")

        # Write out the particles to an HDF5 file

        ctrlCI.timestep_count = 0
        ctrlCI.time = 0.0

        # Run identifier
        ctrlCI.title = "test_ParticleGeneration"
        # Run author
        ctrlCI.author = "tph"

        ### Select output for particles
        ctrlCI.particle_output_file = "particleInitialization.h5part"
        ctrlCI.particle_output_interval = 1
        ctrlCI.particle_output_attributes = ('species_index', 'x', 'y', 'z', 'ux', 'uy', 'uz', 'weight')

        # Check these values
        pCI.check_particle_output_parameters(ctrlCI)

        # Dump the particle data to a file
        pCI.initialize_particle_output_file(ctrlCI)

        # Write the particle attributes
        pCI.write_particle_attributes(ctrlCI)

        return
#    def test_functional_particles(self):ENDDEF


#class TestParticleInitialization(unittest.TestCase):
    def test_something(self):
        """ Check the number of particles in each species.
        """
        fncName = sys._getframe().f_code.co_name
        print '\ntest:', fncName, '('+__file__+')'

        return
#    def test_something(self):ENDDEF

if __name__ == '__main__':
    unittest.main()
