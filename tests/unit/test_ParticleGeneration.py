#!/usr/bin/env python

__version__ = 0.1
__author__ = 'Copyright (C) 2017 L. D. Hughes'
#__all__ = []

import sys
import os
import numpy as np_M
import importlib as im_M
import unittest

import dolfin as df_M

from DT_Module import DTcontrol_C

from Dolfin_Module import *
from UserMesh_FE_XYZ_Module import *

from SegmentedArrayPair_Module import SegmentedArray_C
from Particle_Module import *
from Trajectory_Module import *

from UserUnits_Module import MyPlasmaUnits_C

#STARTCLASS
class TestParticleGeneration(unittest.TestCase):
    """Test user-specified boundary conditions on particles.
    """
    
    def setUp(self):

        # Initialization code common to the tests go here...

        # Common particle inputs
        self.pinCI = ParticleInput_C()

        self.pinCI.precision = np_M.float64
        self.pinCI.particle_integration_loop = 'loop-on-particles'
        self.pinCI.force_precision = np_M.float64

        return

#class TestParticleGeneration:
    def test_1D_particle_source_region(self):
        """Set up a one-dimensional source region and display it.

        """

        ## Set control variables

        ctrlCI = DTcontrol_C()

        ctrlCI.dt = 1.0e-6
        ctrlCI.n_timesteps = 1

        pinCI = self.pinCI

        pinCI.position_coordinates = ['x',] # determines the particle-storage dimensions
        pinCI.force_components = ['x',]

        ### Particle input

        ## Define electron species
        charge = -1.0*MyPlasmaUnits_C.elem_charge
        mass = 1.0*MyPlasmaUnits_C.electron_mass
        dynamics = 'explicit'
        electronsCI = ParticleSpecies_C('electrons', charge, mass, dynamics)

        # Add the electrons to particle input
        pinCI.particle_species = (electronsCI,
                                 )

        ## Make the particle storage array for all species.
        particleCI = Particle_C(pinCI, print_flag=True)

        ## Give the name of the .py file containing special particle data (lists of
        # particles, boundary conditions, source regions, etc.)
        userParticleModule = "UserParticles_1D"

        # Import this module
        UPrt_M = im_M.import_module(userParticleModule)

        particleCI.user_particle_module = userParticleModule
        particleCI.user_particle_class = userParticleClass = UPrt_M.UserParticleDistributions_C
        

        ### Mesh input for a 1D mesh

        # Specify the mesh parameters
        mi1DCI = UserMeshInput_C()
        mi1DCI.pmin = df_M.Point(-10.0)
        mi1DCI.pmax = df_M.Point(10.0)
        mi1DCI.cells_on_side = (20)

        ### Create the 1D particle mesh and add to the Particle_C object
        pmesh1DCI = UserMesh_C(mi1DCI, compute_dictionaries=True, compute_tree=True, plot_flag=False)
        particleCI.pmeshCI = pmesh1DCI

        ### Input for particle sources
        #   For each source:
        #     a. Define the source region
        #     b. Provide a particle-generating function
        #     c. Provide the species name and physical source parameters


        # Create numpy storage needed below
        driftVelocity = np_M.empty(particleCI.particle_dimension, dtype=particleCI.precision)
        ## 1. Hot electron source on left

        # a. Provide the species name and physical parameters
        #    for this source
        speciesName = 'electrons'

        # Check that this species has been defined
        if speciesName in particleCI.species_names:
            charge = particleCI.charge[speciesName]
            mass = particleCI.mass[speciesName]
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
        thermalSpeed = np_M.sqrt(2.0*temp_joule/mass)

        # Set a drift velocity
        driftVelocity[:] = (1.0,)

        # Set desired particles-per-cell
        numberPerCell = 1

        # Collect the parameters for this source in a dictionary
        hotElectronParams = {'species_name': speciesName,
                             'source_distribution_type': sourceDistributionType,
                             'number_density': numberDensity,
                             'thermal_speed': thermalSpeed,
                             'timestep_interval': timeStepInterval,
                             'drift_velocity': driftVelocity,
                             'number_per_cell': numberPerCell}

        # b. Name the particle-creation function.
        maxwellianGenerator = particleCI.add_maxwellian_particles

        # c. Source-region geometry
        xmin = -9.0
        xmax = -4.0
        hotElectronsRegion = RectangularRegion_C(pmesh1DCI, xmin, xmax)

        ## 2. Background electrons over the whole mesh

        ## Specify one or more species to be generated in the 2nd region.

        # a. Provide the species name and physical parameters for this source
        speciesName = 'electrons'

        # Check that this species has been defined above
        if speciesName in particleCI.species_names:
            charge = particleCI.charge[speciesName]
            mass = particleCI.mass[speciesName]
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
        thermalSpeed = np_M.sqrt(2.0*temp_joule/mass)

        # Set a drift velocity
        driftVelocity[:] = (0.5,)

        # Set desired particles-per-cell
        numberPerCell = 10

        # Collect the parameters for this source in a dictionary
        backgroundElectronParams = {'species_name': speciesName,
                                    'source_distribution_type': sourceDistributionType,
                                    'number_density': numberDensity,
                                    'thermal_speed': thermalSpeed,
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

        ### Create input for a particle trajectory object

        # Use an input object to collect initialization data for the trajectory object
        self.trajinCI = TrajectoryInput_C()

        self.trajinCI.maxpoints = None # Set to None to get every point
        self.trajinCI.explicit_dict = {'names': ['x', 'ux',], 'formats': [np_M.float32]*2}

        ## Create the trajectory object and attach it to the particle object.
        # No trajectory storage is created until particles
        # with TRAJECTORY_FLAG on are encountered.
        pCI = particleCI
        trajCI = Trajectory_C(self.trajinCI, ctrlCI, pCI.explicit_species, pCI.implicit_species, pCI.neutral_species)
        particleCI.trajCI = trajCI

        # Add the sources to the particle mesh? or leave them standalone?
#        pmesh1DCI.particle_source_dict = particleSourceDict

        particleCI.particle_source_dict = particleSourceDict

        ## Invoke the source functions and plot the particles
        # time = 0.0
        # for src in particleSourceDict:
        #     print 'Source:', src
        #     (srcRegion, srcFunc, srcParams) = particleSourceDict[src]
        #     srcFunc(time, srcRegion, srcParams)

        
        ctrlCI.timestep_count = 0
        ctrlCI.time = 0.0

        # Run identifier
        ctrlCI.title = "test_ParticleGeneration"
        # Run author
        ctrlCI.author = "tph"

        ### Select output for particles
        ctrlCI.particle_output_file = "p.h5part"
        ctrlCI.particle_output_interval = 1
        ctrlCI.particle_output_attributes = ('species_index', 'x', 'ux', 'weight')

        # Check these values
        particleCI.check_particle_output_parameters(ctrlCI)

        particleCI.add_more_particles(ctrlCI, print_flag=True)        

        # Dump the particle data to a file
        particleCI.initialize_particle_output_file(ctrlCI)

        # Write the particle attributes
        particleCI.write_particle_attributes(ctrlCI)

        return
#    def test_1D_particle_source_region(self):ENDDEF

    def test_2D_particle_source_region(self):
        """Set up a two-dimensional source region and display it.

        """

        from UserMesh_FE_XYZ_Module import UserMeshInput_C

        mi2DCI = UserMeshInput_C()

        mi2DCI.pmin = df_M.Point(-10.0, -10.0)
        mi2DCI.pmax = df_M.Point(10.0, 10.0)
        mi2DCI.cells_on_side = (4, 2)

        ## Boundary conditions for the particles on this mesh

        # Create a 2D particle mesh
#        pmesh2DCI = UserMesh_C(mi2DCI, computeDictionaries=True, computeTree=True, plotFlag=plotFlag)

        return

#class TestParticleGeneration(unittest.TestCase):ENDCLASS

if __name__ == '__main__':
    unittest.main()
