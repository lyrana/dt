#!/usr/bin/env python

# 1D spherically-symmetric expansion of ions and electrons 

__version__ = 0.1
__author__ = 'Copyright (C) 2017 L. D. Hughes'
__all__ = ['Example_C',]

"""
   A spherical source region creates ion-electron pairs.  The
   electrons are hot and the ions are cold.  The electrons pull the
   ions radially outward.
"""

############################## IM. Imported Modules ##############################
import sys
import os
import numpy as np_M
import importlib as im_M
import unittest

import dolfin as df_M

from DT_Module import DTcontrol_C

from Dolfin_Module import *
from UserMesh_y_Fields_Spherical1D_Module import *

from Particle_Module import *
from Trajectory_Module import *

from UserUnits_Module import MyPlasmaUnits_C


############################## SC. Simulation Control ##############################
ctrlObj = DTcontrol_C()
# Run identifier
ctrlObj.title = "sphere1D"
# Run author
ctrlObj.author = "tph"

# Timesteps
ctrlObj.dt = 1.0e-6
ctrlObj.n_timesteps = 1

# Initialize time counters
ctrlObj.timestep_count = 0
ctrlObj.time = 0.0


############################## FM. Field Mesh ##############################

########## FM.1. Provide parameters for a 1D mesh in the radial coordinate ##########
miObj = UserMeshInput_C()

miObj.rmin, miObj.rmax = 1.0, 5.0 # Mesh goes from rmin to rmax in radius
miObj.nr = 10 # Number of divisions in r direction
miObj.stretch = 1.3 # Stretch parameter

########## FM.2. Mark boundaries for BCs on fields ##########

# Name the Dirichlet boundaries and assign integers to them.
# This is done before the mesh is created because sometimes it's
# easier to set boundary conditions first (e.g., if the mesh is
# morphed from a rectangle to a circle)

# Here are the boundary-name -> int pairs used to mark mesh
# facets:
rminIndx = 1
rmaxIndx = 2
fieldBoundaryDict = {'rmin': rminIndx,
                     'rmax': rmaxIndx,
                     }

miObj.field_boundary_dict = fieldBoundaryDict

##FM1 Create the mesh after the particle boundary-conditions are set
## below. See PS.2.1 below.


############################## PS. Particle Species ##############################

# Create an instance of the DTparticleInput class
pinObj = ParticleInput_C()
# Initialize particles
pinObj.precision = np_M.float64
pinObj.particle_integration_loop = 'loop-on-particles'
pinObj.position_coordinates = ['x',] # determines the particle-storage dimensions
pinObj.force_components = ['x',]
pinObj.force_precision = np_M.float64

############### PS.1. Define the species and provide storage ###############

##### PS.1.1. Basic properties and storage #####
speciesName = 'electrons'
charge = -1.0*MyPlasmaUnits_C.elem_charge
mass = 1.0*MyPlasmaUnits_C.electron_mass
dynamics = 'explicit'
electronObj = ParticleSpecies_C(speciesName, charge, mass, dynamics)

speciesName = 'Hplus'
charge = 1.0*MyPlasmaUnits_C.elem_charge
mass = 1.0*MyPlasmaUnits_C.AMU
dynamics = 'explicit'
HplusObj = ParticleSpecies_C(speciesName, charge, mass, dynamics)

# Add the species to particle input:
pinObj.particle_species = (electronObj, HplusObj,
                         )

##### PS.1.2. Make the particle storage array for all species #####
particleObj = Particle_C(pinObj, print_flag=True)

##### PS.1.3. The "UserParticles" module #####
# Give the name of the Python module file containing special particle
# data (lists of particles, boundary conditions, source regions, etc.)
userParticleModule = "UserParticles_1D"

# Import the module
UPrt_M = im_M.import_module(userParticleModule)

# Add it to the particle object
particleObj.user_particle_module = userParticleModule
particleObj.user_particle_class = userParticleClass = UPrt_M.UserParticleDistributions_C

############### PS.2. Particle boundary conditions, mesh, sources ###############

##### PS.2.1. Mark boundaries to apply BCs on particles #####
# In some cases, it's easier to mark boundaries before the final mesh is created.

# These are the (boundary-name, int) pairs used to mark mesh facets:
rminIndx = 1
rmaxIndx = 2
particleBoundaryDict = {'rmin': rminIndx,
                        'rmax': rmaxIndx,
                        }

# Add the above dict to the mesh input
miObj.particle_boundary_dict = particleBoundaryDict

##### PS.2.2. Make the mesh and add it to the particle object #####

# First complete the mesh setup now that particle boundary-conditions have been
# set (see #FM1 above).
plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name
meshObj = UserMesh_C(meshInputCI=miObj, compute_dictionaries=True, compute_tree=False, plot_flag=True, plot_title=plotTitle)

# In this simulation, the particle mesh is the same object as the field mesh
particleObj.pmeshObj = meshObj

########## PS.3. Particle source regions ##########

        ### Input for particle sources
        #   For each source:
        #     a. Define the source region
        #     b. Provide a particle-generating function
        #     c. Provide the species name and physical source parameters


##### PS.3.1. Electron source at r=0 #####
speciesName = 'electrons'

# Check that this species has been defined
if speciesName in particleObj.species_names:
    charge = particleObj.charge[speciesName]
    mass = particleObj.mass[speciesName]
else:
    print "The species", speciesName, "has not been defined"
    sys.exit()

## PS.3.1.1. Define the electron distribution function & parameters ##

# Name the particle-creation function:
particleGenerator = particleObj.create_maxwellian_particles
# !!!The function could provide it's own type
#sourceDistributionType = 'functional'

# Provide the charge-density creation rate:
chargeDensityRate = -1.0

# Compute a number-density rate from the charge-density rate
timeStepInterval = 1 # Timesteps between invocations
numberDensity = chargeDensityRate*timeStepInterval*ctrlObj.dt/charge
# Check for positivity
if numberDensity < 0:
    print "Check number density for species", speciesName, "is negative. Should be positive"
    sys.exit()

# Compute a value for thermalSpeed from the temperature in eV
temperature_eV = 2.0

temp_joule = temperature_eV*MyPlasmaUnits_C.elem_charge
thermalSpeed = np_M.sqrt(2.0*temp_joule/mass)

# Set a drift velocity
vdrift_x = 1.0
vdrift_y = 2.0
vdrift_z = 3.0

driftVelocity = (vdrift_x, vdrift_y, vdrift_z)

# Set desired particles-per-cell
numberPerCell = 1

## PS.3.1.2. Collect the parameters for this source in a dictionary
electronSourceParams = {'species_name': speciesName,
#                     'source_distribution_type': sourceDistributionType,
                     'number_density': numberDensity,
                     'thermal_speed': thermalSpeed,
                     'timestep_interval': timeStepInterval,
                     'drift_velocity': driftVelocity,
                     'number_per_cell': numberPerCell}


## PS.3.1.3. Source geometry for electrons
xmin = -9.0
xmax = -4.0
electronSourceRegion = RectangularRegion_C(particleObj.pmeshObj, xmin, xmax)

##### PS.3.2. Hplus source at r=0 #####
speciesName = 'Hplus'

# Check that this species has been defined
if speciesName in particleObj.species_names:
    charge = particleObj.charge[speciesName]
    mass = particleObj.mass[speciesName]
else:
    print "The species", speciesName, "has not been defined"
    sys.exit()

## PS.3.2.1. Define the Hplus distribution function & parameters ##

# Name the particle-creation function:
particleGenerator = particleObj.create_maxwellian_particles
# !!!The function could provide it's own type
#sourceDistributionType = 'functional'

# Provide the charge-density creation rate:
chargeDensityRate = 1.0

# Compute a number-density rate from the charge-density rate
timeStepInterval = 1 # Timesteps between invocations
numberDensity = chargeDensityRate*timeStepInterval*ctrlObj.dt/charge
# Check for positivity
if numberDensity < 0:
    print "Check number density for species", speciesName, "is negative. Should be positive"
    sys.exit()

# Compute a value for thermalSpeed from the temperature in eV
temperature_eV = 0.1 # About 1000K

temp_joule = temperature_eV*MyPlasmaUnits_C.elem_charge
thermalSpeed = np_M.sqrt(2.0*temp_joule/mass)

# Set a drift velocity
vdrift_x = 0.0
vdrift_y = 0.0
vdrift_z = 0.0

driftVelocity = (vdrift_x, vdrift_y, vdrift_z)

# Set desired particles-per-cell
numberPerCell = 1

## PS.3.2.2. Collect the parameters for this source in a dictionary ##
HplusSourceParams = {'species_name': speciesName,
#                     'source_distribution_type': sourceDistributionType,
                     'number_density': numberDensity,
                     'thermal_speed': thermalSpeed,
                     'timestep_interval': timeStepInterval,
                     'drift_velocity': driftVelocity,
                     'number_per_cell': numberPerCell}


## PS.3.2.3. Source geometry for Hplus ##
xmin = -9.0
xmax = -4.0
HplusSourceRegion = RectangularRegion_C(particleObj.pmeshObj, xmin, xmax)


##### PS.3.3. Put all the particle source data into a dictionary #####

# The dictionary keys are mnemonic source names.
# The dictionary values are lists of tuples giving the 
#   a. the physical source parameters.
#   b. particle creation function and
#   c. source region

## Note: the names here are source-region names, NOT species names!
particleSourceDict = {'electron_source': (electronSourceParams, particleGenerator, electronSourceRegion),
                      'Hplus_source': (HplusSourceParams, particleGenerator, HplusSourceRegion),
                      }

particleObj.particle_source_dict = particleSourceDict


############################## FS. Field Solver ##############################

########## FS.1. The scalar electric potential ##########

##### FS.1.1. Boundary values of the potential #####
phiVals = {'rmin':  0.0,
           'rmax': -1.5,
           }

## Save these in a boundary-condition dictionary
phiBCs = dict( (bnd, [fieldBoundaryDict[bnd], phiVals[bnd]]) for bnd in fieldBoundaryDict.keys())

##### FS.1.2. Properties of the electric potential trial function #####
phiElementType = 'Lagrange'
phiElementDegree = 2
phiFieldType = 'scalar'

##### FS.1.3. Storage for the electric potential solution
phi = df_M.Field_C(meshObj=meshObj,
                   element_type=phiElementType,
                   element_degree=phiElementDegree,
                   field_type=phiFieldType)

########## FS.2. The electric field (gradient of scalar potential) ##########

##### FS.2.1. Properties of the electric field #####
if phiElementDegree == 1:
    # For linear elements, grad(phi) is discontinuous across
    # elements. To represent this field, we need Discontinuous Galerkin
    # elements.
    electricFieldElementType = 'DG'
else:
    electricFieldElementType = 'Lagrange'

##### FS.2.2. Storage for the electric field #####
negElectricField = df_M.Field_C(meshObj=meshObj,
                              element_type=electricFieldElementType,
                              element_degree=phiElementDegree-1,
                              field_type='vector')

########## FS.3.1. The Poisson solver parameters ##########
linearSolver = 'lu'
preconditioner = None

########## FS.3.2. Create the Poisson solver object ##########
fieldBoundaryMarker = meshObj.field_boundary_marker
poissonsolveObj = UserPoissonSolve_C(phi,
                                   linearSolver, preconditioner,
                                   fieldBoundaryMarker, phiBCs,
                                   neg_electric_field=negElectricField)

########## FS.4. The source terms for the electric field ##########

##### FS.4.1 Charge-density obtained from kinetic particles #####

## FS.4.1.1 Create number-density arrays on the particle mesh for each kinetic species ##
number_density_element_type = 'Lagrange'
number_density_element_degree = 1
number_density_field_type = 'scalar'

for s in particleObj.species_names:
    particleObj.number_density_dict[s] = Field_C(meshCI=particleObj.pmeshObj,
                                               element_type=number_density_element_type,
                                               element_degree=number_density_element_degree,
                                               field_type=number_density_field_type)

## FS.4.1.2 Create a charge-density array on the field mesh ##
chargeDensity = Field_C(meshCI=meshObj,
                         element_type=number_density_element_type,
                         element_degree=number_density_element_degree,
                         field_type=number_density_field_type)

# Attach this to the particle object?
#particleObj.charge_density = chargeDensity

############################## INIT. Initialization ##############################

########## INIT.1. Do an initial field solve without particles ##########

poissonsolveObj.solve_for_phi(plot_flag=plotFlag, plot_title=plotTitle)

########## INIT.2. Initialize the particles ##########

##### INIT.2.1. Loop on initialized particle species
if particleObj.initial_particles_dict is not None:
    printFlags = {}
    for s in particleObj.species_names: printFlags[s] = False
    particleObj.initialize_particles(printFlags)

    ##### INIT.2.2. Recompute the initial electrostatic field in the presence of these particles #####
    particleObj.accumulate_charge_density_from_particles(chargeDensity)
    poissonsolveObj.solve_for_phi(plot_flag=plotFlag, plot_title=plotTitle)
else:
    print "There are no initial particles in this simulation"

          
############################## RUN. Integrate forward in time ##############################

for istep in xrange(ctrlObj.n_timesteps):

    ########## RUN.1. Advance the particles one timestep ##########
    particleObj.move_particles_in_electrostatic_field(ctrlObj.dt, neg_electric_field)

    ########## RUN.2. Add particles from source regions ##########
    if particleObj.particle_source_dict is not None:
        particleObj.add_more_particles(ctrlObj, print_flag=True)        

    ########## RUN.3. Update the electric charge density ##########
    particleObj.compute_charge_density(chargeDensity)

    ########## RUN.4. Update the electric field ##########
    plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name
    poissonsolveObj.solve_for_phi(plot_flag=plotFlag, plot_title=plotTitle)


    ctrlObj.time_step += 1
    ctrlObj.time += ctrlObj.dt
