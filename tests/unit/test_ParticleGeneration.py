#!/usr/bin/env python

__version__ = 0.1
__author__ = 'Copyright (C) 2017 L. D. Hughes'
#__all__ = []

import sys
import os
import numpy as np_m
import importlib as im_m
import unittest

import dolfin as df_m

from DT_Module import DTcontrol_C

from Dolfin_Module import *
from UserMesh_y_Fields_FE_XYZ_Module import *

from SegmentedArrayPair_Module import SegmentedArray_C
from Particle_Module import *
from Trajectory_Module import *

from UserUnits_Module import MyPlasmaUnits_C

#STARTCLASS
class TestParticleGeneration(unittest.TestCase):
    """Test particle-generation functions.
    """
    
    def setUp(self):

        # Initialization code common to the tests go here...

        # Common particle inputs
        self.pin = ParticleInput_C()

        self.pin.precision = np_m.float64
        self.pin.particle_integration_loop = 'loop-on-particles'
        self.pin.force_precision = np_m.float64

        return

#class TestParticleGeneration:
    def test_1D_particle_source_region(self):
        """Set up a one-dimensional source region and display it.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest: ', fncName, '('+__file__+')'

        ## Set control variables

        ctrlCI = DTcontrol_C()

        ctrlCI.dt = 1.0e-6
        ctrlCI.n_timesteps = 1

        pin = self.pin

        pin.position_coordinates = ['x',] # determines the particle-storage dimensions
        pin.force_components = ['x',]

        ### Particle input

        ## Define electron species
        charge = -1.0*MyPlasmaUnits_C.elem_charge
        mass = 1.0*MyPlasmaUnits_C.electron_mass
        dynamics = 'explicit'
        electronsCI = ParticleSpecies_C('electrons', charge, mass, dynamics)

        # Add the electrons to particle input
        pin.particle_species = (electronsCI,
                                 )

        ## Make the particle storage array for all species.
        particle_P = Particle_C(pin, print_flag=True)

        ## Give the name of the .py file containing special particle data (lists of
        # particles, boundary conditions, source regions, etc.)
        userParticleModule = "UserParticles_1D"

        # Import this module
        UPrt_M = im_m.import_module(userParticleModule)

        particle_P.user_particle_module = userParticleModule
        particle_P.user_particle_class = userParticleClass = UPrt_M.UserParticleDistributions_C
        

        ### Mesh input for a 1D mesh

        # Specify the mesh parameters
        umi1DCI = UserMeshInput_C()
        umi1DCI.pmin = df_m.Point(-10.0)
        umi1DCI.pmax = df_m.Point(10.0)
        umi1DCI.cells_on_side = (20,) # Need the comma to indicate a tuple

        ### Create the 1D particle mesh and add to the Particle_C object
        pmesh1DCI = UserMesh_C(umi1DCI, compute_dictionaries=True, compute_tree=True, plot_flag=False)
        particle_P.pmeshCI = pmesh1DCI

        ### Input for particle sources
        #   For each source:
        #     a. Define the source region
        #     b. Provide a particle-generating function
        #     c. Provide the species name and physical source parameters


        # Create numpy storage needed below
#        driftVelocity = np_m.empty(particle_P.particle_dimension, dtype=particle_P.precision)
        ## 1. Hot electron source on left side of the mesh

        # a. Provide the species name and physical parameters
        #    for this source
        speciesName = 'electrons'

        # Check that this species has been defined
        if speciesName in particle_P.species_names:
            charge = particle_P.charge[speciesName]
            mass = particle_P.mass[speciesName]
        else:
            print "The species", speciesName, "has not been defined"
            sys.exit()

        sourceDistributionType = 'functional'
        # Compute a physical number-density creation rate for this source region
        # The value should always be positive 
        chargeDensityRate = -1.0
        # The timestep interval between calls to the creation function
        timeStepInterval = 1
        # Get number-density added per invocation
        numberDensity = chargeDensityRate*timeStepInterval*ctrlCI.dt/charge
        # Check for positivity
        if numberDensity < 0:
            print "Check number density for species", speciesName, "is negative. Should be positive"
            sys.exit()

        # Compute a value for thermalSpeed from the temperature in eV
        temperature_eV = 2.0
        temp_joule = temperature_eV*MyPlasmaUnits_C.elem_charge
        thermalSpeed = np_m.sqrt(2.0*temp_joule/mass)

        velocity_coordinate_system = 'x_y_z' # This sets the interpretation of the
                                             # following values
        # Set a drift velocity
        vdrift_1 = 1.0
        vdrift_2 = 2.0
        vdrift_3 = 3.0

        driftVelocity = (vdrift_1, vdrift_2, vdrift_3)

        # Set desired particles-per-cell
        numberPerCell = 1

        # Collect the parameters for this source in a dictionary
        hotElectronParams = {'species_name': speciesName,
                             'source_distribution_type': sourceDistributionType,
                             'number_density': numberDensity,
                             'thermal_speed': thermalSpeed,
                             'velocity_coordinate_system': velocity_coordinate_system,
                             'drift_velocity': driftVelocity,
                             'timestep_interval': timeStepInterval,
                             'number_per_cell': numberPerCell}

        # b. Name the particle-creation function.
        maxwellianGenerator = particle_P.create_maxwellian_particles

        # c. Source-region geometry
        xmin = -9.0
        xmax = -4.0
        hotElectronsRegion = RectangularRegion_C(pmesh1DCI, xmin, xmax)

        ## 2. Background electrons over the whole mesh

        ## Specify one or more species to be generated in the 2nd region.

        # a. Provide the species name and physical parameters for this source
        speciesName = 'electrons'

        # Check that this species has been defined above
        if speciesName in particle_P.species_names:
            charge = particle_P.charge[speciesName]
            mass = particle_P.mass[speciesName]
        else:
            print "The species", speciesName, "has not been defined"
            sys.exit()

        sourceDistributionType = 'functional'
        # Compute a value for numberDensity, the physical number-density
        # creation rate for this source region.
        # The value should always be positive.
        chargeDensityRate = -0.1
        # The timestep interval between calls to the creation function
        timeStepInterval = 1
        # Get charge-density per invocation
        numberDensity = chargeDensityRate*timeStepInterval*ctrlCI.dt/charge
        # Check for positivity
        if numberDensity < 0:
            print "Check number density for species", speciesName, "is negative. Should be positive"
            sys.exit()

        # Compute a value for thermalSpeed
        temperature_eV = 2.0
        temp_joule = temperature_eV*MyPlasmaUnits_C.elem_charge
        thermalSpeed = np_m.sqrt(2.0*temp_joule/mass)

        velocity_coordinate_system = 'x_y_z' # This sets the interpretation of the
                                             # following values
        # Set a drift velocity
        vdrift_1 = 0.5
        vdrift_2 = 1.5
        vdrift_3 = 2.5

        driftVelocity = (vdrift_1, vdrift_2, vdrift_3)

        # Set desired particles-per-cell
        numberPerCell = 10

        # Collect the parameters for this source in a dictionary
        backgroundElectronParams = {'species_name': speciesName,
                                    'source_distribution_type': sourceDistributionType,
                                    'number_density': numberDensity,
                                    'thermal_speed': thermalSpeed,
                                    'velocity_coordinate_system': velocity_coordinate_system,
                                    'drift_velocity': driftVelocity,
                                    'timestep_interval': timeStepInterval,
                                    'number_per_cell': numberPerCell}

        # b. Name the first function that will generate particles in the above region:
        #      Just use same maxwellianGenerator function as above.

        # c. Source geometry: in this case, it's the whole mesh:
        wholeMesh = WholeMesh_C(pmesh1DCI)

        ## Put all the particle source data into a dictionary

        # The dictionary keys are mnemonic source names.
        # The dictionary values are lists of tuples giving the 
        #   a. the physical source parameters.
        #   b. particle creation function and
        #   c. source region

        ## Note: the names here are source-region names, NOT species names!
        particleSourceDict = {'hot_electron_source': (hotElectronParams, maxwellianGenerator, hotElectronsRegion),
                              'background_electron_source': (backgroundElectronParams, maxwellianGenerator, wholeMesh),
                              }

        particle_P.particle_source_dict = particleSourceDict

        ### Create input for a particle trajectory object

        # Use an input object to collect initialization data for the trajectory object
        self.trajinCI = TrajectoryInput_C()

        self.trajinCI.maxpoints = None # Set to None to get every point
        self.trajinCI.explicit_dict = {'names': ['x', 'ux',], 'formats': [np_m.float32]*2}

        ## Create the trajectory object and attach it to the particle object.
        # No trajectory storage is created until particles
        # with TRAJECTORY_FLAG on are encountered.
        p_P = particle_P
        trajCI = Trajectory_C(self.trajinCI, ctrlCI, p_P.explicit_species, p_P.implicit_species, p_P.neutral_species)
        particle_P.trajCI = trajCI

        ## Invoke the source functions and write out the particles
        
        ctrlCI.timeloop_count = 0
        ctrlCI.time = 0.0

        # Run identifier
        ctrlCI.title = "test_ParticleGeneration"
        # Run author
        ctrlCI.author = "tph"

        ### Select output for particles
        ctrlCI.particle_output_file = "particleGeneration.h5part"
        ctrlCI.particle_output_interval = 1
        ctrlCI.particle_output_attributes = ('species_index', 'x', 'ux', 'weight')

        # Check these values
        particle_P.check_particle_output_parameters(ctrlCI)

        particle_P.add_more_particles(ctrlCI, print_flag=True)        

        # Dump the particle data to a file
        particle_P.initialize_particle_output_file(ctrlCI)

        # Write the particle attributes
        particle_P.write_particles_to_file(ctrlCI)

        return
#    def test_1D_particle_source_region(self):ENDDEF

    def test_2D_particle_source_region(self):
        """Set up a two-dimensional source region and display it.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest: ', fncName, '('+__file__+')'

#        from UserMesh_y_Fields_FE_XYZ_Module import UserMeshInput_C

        umi2DCI = UserMeshInput_C()

# Change these to TUPLES instead of Points (avoid df_m at toplevel)
        umi2DCI.pmin = df_m.Point(-10.0, -10.0)
        umi2DCI.pmax = df_m.Point(10.0, 10.0)
        umi2DCI.cells_on_side = (4, 2)

        ## Boundary conditions for the particles on this mesh

        # Create a 2D particle mesh
#        pmesh2DCI = UserMesh_C(umi2DCI, computeDictionaries=True, computeTree=True, plotFlag=plotFlag)

        return
#    def test_2D_particle_source_region(self):ENDDEF

#class TestParticleGeneration(unittest.TestCase):ENDCLASS

if __name__ == '__main__':
    unittest.main()
