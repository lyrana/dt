#!/usr/bin/env python

# 1D spherically-symmetric expansion of ions and electrons
# See sphere1D.ods for setting up parameters

__version__ = 0.1
__author__ = 'Copyright (C) 2017 L. D. Hughes'
__all__ = ['Example_C',]

"""A spherical source region creates ion-electron pairs.  The electrons are hot and
   the ions are cold, but may be given a radial velocity.  The electrons pull the
   ions radially outward.

"""
#### Commentary
#
#  The key parameters are:
#
#      1. The electron temperature ('T_e_eV'). This determines how fast the electrons
#         escape.
#
#      2. The potential difference ('confining_potential') between the center and the
#         edge of the region.  This is the confining potential for the electrons.
#      
#

#### Section Index

### To search for a section: type a space after the last period of the section
### number below, e.g., PS.<space>, or PS.1.3.<space>

## IM. Imported Modules
## OF. Terminal output flags
## PC. Physical constants
## SP. Scratch-pad for calculating parameter values.
##     SP.1. Calculate spatial parameters
##     SP.2. Calculate electron parameters
##     SP.3. Calculate ion parameters
##     SP.4. Calculate field parameters
## SC. Simulation Control
##     Set timestep value and number of timesteps.
##     Output: Specify field and particle output files and output frequency.
## FM. Field Mesh
## PS. Particle Species
##     PS.1. Define the species and provide storage.
##           PS.1.2. Allocate particle storage.
##           PS.1.3. Specify the name of the UserParticlesModule (Deprecate)
##     PS.2. Particle boundary conditions, particle mesh
##           PS.2.1. Mark boundaries to apply BCs on particles
##           PS.2.2. Make the particle mesh and add it to the particle object
##           PS.2.3. Connect the boundary-condition callback functions
##     PS.3. Particle source regions
##           PS.3.1. Hplus source near r=0
##           PS.3.2. Electron source near r=0
##           PS.3.3. Put all the particle source data into a dictionary
## TRJ. Particle Trajectories
## SRC. Source terms for the electric potential
## FS. Field Solver
##     FS.1. The scalar electric potential
##           FS.1.1. Set boundary values of the potential
##     FS.2. The electric field (gradient of scalar potential)
##           FS.2.1. Properties of the electric field
##           FS.2.2. Storage for the electric field
##     FS.3. The potential field solver
##           FS.3.1. The potential field solver algorithm parameters
##           FS.3.2. Create the potential field solver object
## HST. Time-Histories
## CHK. Check the input
## INIT. Initialization
##       INIT.1. Do an initial field solve without particles
##       INIT.2. Recompute initial fields including initial particle spacecharge
##               and random electric field
##               INIT.2.1. Recalculate initial electric field
##                         INIT.2.1.1. Recompute the initial electrostatic field
##                                     with particle spacecharge
##                         INIT.2.1.2. Add in a random electric field
##               INIT.2.2. Record initial history data
## RUN. Integrate forward in time
##      RUN.0. Increment the time counters
##      RUN.1. Advance the particles one timestep
##      RUN.2. Add particles from source regions
##      RUN.3. Write a snapshot of the particles
##      RUN.4. Update the electric charge density
##      RUN.5. Update the electric field
##             RUN.5.1. Set field-plotting parameters
##             RUN.5.2. Call the field solver
##             RUN.5.3. Add an external electric field
##             RUN.5.4. Write a snapshot of the fields
##      RUN.6. Record history data for this timestep
## FIN. Finish the simulation
##      FIN.1. Write the last particle data and close the output file.
##      FIN.2. Record the LAST point on the particle trajectories.
##      FIN.3. Plot the particle trajectories on the particle mesh.
##      FIN.4. Plot the particle trajectories in phase-space.
##      FIN.5. Plot the recorded history data.
##      FIN.6. Write final log message


############################## IM. Imported Modules ##############################
import sys
import os
import numpy as np_m
import importlib as im_M
import unittest

import dolfin as df_m

from DT_Module import DTcontrol_C, DTscratch_C
from Dolfin_Module import *
from UserMesh_y_Fields_Spherical1D_Module import * # User input for mesh and fields

from UserParticles_1D import * # User input for particles
from Particle_Module import *
from RecordedData_Module import *

from UserUnits_Module import MyPlasmaUnits_C

fileName = __file__+':'

############################## OF. Terminal output flags ##############################

emitInput = True # Write input values to standard output.
emitInputOnly = False # Write input values to standard output and then quit.
pauseAfterEmit = True # Wait for user input before continuing.  Can be changed to
                      # False by using ctrl.get_keyboard_input()

pauseAfterTimestep = True # Wait for user input before continuing to the next
                          # timestep.  Can be changed to False by using
                          # ctrl.get_keyboard_input()
                      
plotInitialFields = True
plotFieldsEveryTimestep = True # Plot phi and -E each timestep.
plotFinalFields = True

plotTrajectoriesOnMesh = True # Plot the trajectory particles on top of the mesh.
plotTrajectoriesInPhaseSpace = True # Plot the 2D phase-space of trajectory particles.

plotTimeHistories = True # Plot the recorded values vs. time.

############################## PC. Physical constants ##############################

m_e = 1.0*MyPlasmaUnits_C.electron_mass
q_e = 1.0*MyPlasmaUnits_C.elem_charge

m_Hplus = 1.0*MyPlasmaUnits_C.AMU

eps_0 = MyPlasmaUnits_C.epsilon_0
k_B = MyPlasmaUnits_C.boltzmann_constant


############################## SP. Scratch-pad ##############################

# Use the scratch-pad to compute values needed to set up the simulation

if emitInput is True:
    print ""
    print "********** Scratch-pad for computing simulation parameters **********"

# Create a class for scratch-pad variables
scr = DTscratch_C()

########## SP.1. Calculate spatial parameters

# Start with the size of the computational domain
#scr.rmin = 1.0
scr.rmin = 0.0
scr.rmax = 4.0
# Set the spatial resolution
scr.nr = 10

# Compute a rough cell-size (The mesh can be stretched non-uniformly)
scr.dr_av = (scr.rmax-scr.rmin)/scr.nr

# Specify the desired Debye length as a fraction of the cell-size
scr.lambda_D = 0.5*scr.dr_av
print "Electron Debye length = %.3g m" % scr.lambda_D
ctrl = DTcontrol_C()
if emitInput is True:
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

########## SP.2. Calculate electron parameters

# Specify the electron temperature
scr.T_e_eV = 6.0

# Convert this to an electron thermal velocity
scr.vthermal_e = np_m.sqrt(q_e*scr.T_e_eV/m_e)
print "Electron thermal velocity = %.3g m/s" % scr.vthermal_e
if emitInput is True:
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

# Compute the electron plasme frequency needed to get the above Debye length
scr.omega_e = scr.vthermal_e/scr.lambda_D
print "Electron plasma frequency = %.3g rad/s" % scr.omega_e
if emitInput is True:
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

# Compute the explicit timestep limit from the electron plasma frequency
scr.dt_max = 2.0/scr.omega_e
print "Maximum stable timestep = %.3g s" % scr.dt_max
if emitInput is True:
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

# Compute the electron density from the plasma frequency
scr.electron_density = m_e*eps_0*(scr.omega_e/q_e)**2
print "Electron density = %.3g per m^{-3}" % scr.electron_density
if emitInput is True:
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

########## SP.3. Calculate ion parameters

# Compute the ion sound speed
scr.ion_charge_state = 1.0
scr.ion_sound_speed = np_m.sqrt(scr.ion_charge_state*q_e*scr.T_e_eV/m_Hplus)
print "Ion sound speed = %.3g m/s" % scr.ion_sound_speed
if emitInput is True:
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

#electron_plasma_frequency = q_e*np_m.sqrt(n_e/(m_e*eps_0))

########## SP.4. Set and calculate field parameters

# Set the confining potential (in Volts)
scr.confining_potential = 10.0*scr.T_e_eV
scr.spacecharge_factor = 0.01

scr.external_electric_field_amplitude = 0.1

############################## SC. Simulation Control ##############################

### Titles
# Run identifier
ctrl.title = "sphere1D"
# Run author
ctrl.author = "tph"

### Timestepping
ctrl.n_timesteps = 100 # 100

ctrl.dt = 4.0e-7 # from sphere1D.ods
if ctrl.dt > scr.dt_max:
    print "%s\n(DnT WARNING) Timestep exceeds stability limit by %.3g" % (fileName, ctrl.dt/scr.dt_max)
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)


# Initialize time counters
ctrl.timeloop_count = 0
ctrl.time = 0.0

# Set random seed
ctrl.random_seed = 1
np_m.random.seed(ctrl.random_seed)

### External electric field switch
ctrl.apply_external_electric_field = True

### Set parameters for writing field output files
ctrl.potential_field_output_file = ctrl.title + "-phi.pvd"
ctrl.electric_field_output_file = ctrl.title + "-E.pvd"
ctrl.field_output_interval = 1

### Set parameters for writing particle output to an h5part file
#ctrl.particle_output_file = "sphere1D.h5part"
ctrl.particle_output_file = ctrl.title + ".h5part"
ctrl.particle_output_interval = 1
# Available attributes: 'species_index' and any particle-record name
ctrl.particle_output_attributes = ('species_index', 'x', 'ux', 'unique_ID')

if emitInput is True:
    print ""
    print "********** Simulation Control Attributes **********"
    print ctrl
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)


#    print ctrl.__dict__

#output = DToutput_C()


############################## FM. Field Mesh ##############################

########## FM.1. Provide parameters for a 1D mesh in the radial coordinate
#from UserMesh_y_Fields_Spherical1D_Module import *
umi = UserMeshInput1DS_C()

umi.rmin, umi.rmax = scr.rmin, scr.rmax # Mesh goes from rmin to rmax in radius
umi.nr = scr.nr # Number of divisions in r direction
umi.stretch = 1.3 # Stretch parameter

#Note: more attributes are added below

########## FM.2. Mark boundaries for BCs on fields

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

umi.field_boundary_dict = fieldBoundaryDict

if emitInput is True:
    print ""
    print "********** Field Boundary Dictionary (fieldBoundaryDict.__dict__) **********"
#    print "fieldBoundaryDict =", fieldBoundaryDict
    for k, v in fieldBoundaryDict.iteritems():
        if type(v) is np_m.float64:
            print " %s = %.3g" % (k, v)
        else:
            print '', k, '=', v
#        print '', k, '=', v
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

##FM1 Create the mesh after the particle boundary-conditions are set
## below. See PS.2.1. below.


############################## PS. Particle Species ##############################

# Create an instance of the DTparticleInput class
pin = ParticleInput_C()
# Initialize particles
pin.precision = np_m.float64
pin.particle_integration_loop = 'loop-on-particles'
pin.coordinate_system = '1D-spherical-radius'
pin.position_coordinates = ['x',] # determines the particle-storage dimensions. This
                                  # is doubled to get the phase-space coordinates
pin.force_components = ['x',]
pin.force_precision = np_m.float64

############### PS.1. Define the species and provide storage ###############

##### PS.1.1. Basic properties and storage
speciesName = 'electron'
charge = -1.0*q_e
mass = 1.0*m_e
dynamics = 'explicit'
electron_S = ParticleSpecies_C(speciesName, charge, mass, dynamics)

speciesName = 'Hplus'
charge = 1.0*q_e
mass = 1.0*m_Hplus
dynamics = 'explicit'
Hplus_S = ParticleSpecies_C(speciesName, charge, mass, dynamics)

# Add the two species to particle input:
pin.particle_species = (electron_S, Hplus_S,)

if emitInput is True:
    print ""
    print "********** Particle Species Attributes (particle_species.__dict__) **********"
    for s in pin.particle_species:
        for k, v in s.__dict__.iteritems():
            if type(v) is np_m.float64:
                print " %s = %.3g" % (k, v)
            else:
                print '', k, '=', v
#            print '', k, '=', v
        print ""
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

if emitInput is True:
    print ""
    print "********** Particle Input Attributes (pin.__dict__) **********"
#    print "Particle Input =", pin.__dict__
    for k, v in pin.__dict__.iteritems():
            if type(v) is np_m.float64:
                print " %s = %.3g" % (k, v)
            else:
                print '', k, '=', v
#        print '', k, '=', v
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

##### PS.1.2. Make particle storage for the species defined above
particle_P = Particle_C(pin, print_flag=True)

##### PS.1.3. Name the particular "UserParticles" module to use
# Give the name of the Python module file containing special particle
# data (lists of particles, boundary conditions, source regions, etc.)
#userParticlesModuleName = "UserParticles_1D"
# See Imported Modules above

# Import the module
#userParticlesModule = im_M.import_module(userParticlesModuleName)

# Add it to the particle object
#?needed? particle_P.user_particles_module = userParticlesModuleName
##particle_P.user_particles_class = userParticlesClass = userParticlesModule.UserParticleDistributions_C

############### PS.2. Particle boundary conditions and mesh ###############

##### PS.2.1. Mark boundaries to apply BCs on particles
# In some cases, it's easier to mark boundaries before the final mesh is created.

# These are the (boundary-name, int) pairs used to mark mesh facets:
rminIndx = 1
rmaxIndx = 2
particleBoundaryDict = {'rmin': rminIndx,
                        'rmax': rmaxIndx,
                        }

if emitInput is True:
    print ""
    print "********** Particle Boundary Dictionary (particleBoundaryDict) **********"
#    print "particleBoundaryDict =", particleBoundaryDict
    for k, v in particleBoundaryDict.iteritems():
        if type(v) is np_m.float64:
            print " %s = %.3g" % (k, v)
        else:
            print '', k, '=', v
#        print '', k, '=', v
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

# Add the above dict to the mesh input
umi.particle_boundary_dict = particleBoundaryDict

if emitInput is True:
    print ""
    print "********** Field Mesh Input Attributes **********"
    print umi
#    print umi.__dict__
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

##### PS.2.2. Make the mesh and add it to the particle object

# First complete the mesh setup now that particle boundary-conditions have been
# set (see #FM1 above).
if emitInput is True:
    print ""
    print "********** Construct and plot the mesh **********"
plotTitle = os.path.basename(__file__) + ": "
mesh_M = UserMesh1DS_C(umi, compute_dictionaries=True, compute_tree=False, plot_flag=emitInput, plot_title=plotTitle)

# In this simulation, the particle mesh is the same object as the field mesh
particle_P.pmesh_M = mesh_M

##### PS.2.3. Connect the boundary-condition callback functions

# See UserParticleBoundaryFunctions_C in userParticlesModule (defined above) for the
# definitions of the facet-crossing callback functions.
#userPBndFns = userParticlesModule.UserParticleBoundaryFunctions_C(particle_P.position_coordinates, particle_P.dx)
userPBndFns = UserParticleBoundaryFunctions_C(particle_P.position_coordinates, particle_P.dx)

# Make the particle-mesh boundary-conditions object and add it to the particle
# object.
spNames = particle_P.species_names
pmeshBCs = ParticleMeshBoundaryConditions_C(spNames, mesh_M, userPBndFns, print_flag=False)
particle_P.pmesh_bcs = pmeshBCs


########## PS.3. Particle source regions

        ### Input for particle sources
        #   For each source 'n':
        #     PS.3.n.1 Give name and physical properties of the species for this source
        #     PS.3.n.2 Give the source-region geometry
        #     PS.3.n.3 Provide the particle-generating function and
        #              the physical parameters
        #     PS.3.n.4 Collect the input for this source into a dictionary
        #   Collect all the sources in a dictionary

##### PS.3.1. Hplus source near r=0

## PS.3.1.1. The species created in this source
ionSpeciesName = 'Hplus'
# Check that this species has been defined
if ionSpeciesName in particle_P.species_names:
    ionCharge = particle_P.charge[ionSpeciesName]
    ionMass = particle_P.mass[ionSpeciesName]
else:
    errorMsg = "(DnT ERROR) %s The species %s has not been defined" % (fncName, ionSpeciesName)
    sys.exit(errorMsg)

## PS.3.1.2. Geometry of Hplus source region
rminHplus = scr.rmin
rmaxHplus = rminHplus + 2.0*scr.dr_av
HplusSourceRegion = RectangularRegion_C(particle_P.pmesh_M, rminHplus, rmaxHplus)

if emitInput is True:
    print ""
    print "********** Hplus Source Region Attributes **********"
#    print "HplusSourceRegion =", HplusSourceRegion.__dict__
    print HplusSourceRegion
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

## PS.3.1.3. Define the Hplus generating function & parameters

# Name the particle-generating function:
particleGenerator = particle_P.create_maxwellian_particles
# !!!The function could provide it's own type
#sourceDistributionType = 'functional'

# Set the drift velocity >= the ion sound speed
velocity_coordinate_system = 'r_theta_z' # This sets the interpretation of the
                                         # following values
# vdrift_r = 2.64e4 # 1.1 times the sound speed
# vdrift_theta = 0.0
# vdrift_z = 0.0
#vdrift_1 = 1.1*scr.ion_sound_speed # 1.1 times the sound speed
vdrift_1 = 10.0*scr.ion_sound_speed
vdrift_2 = 0.0
vdrift_3 = 0.0

# Compute the flux of ions per unit area out of the source region (assuming the ion
# charge-density equals the electron charge-density):
ion_flux = scr.electron_density*vdrift_1/scr.ion_charge_state

# The ions are replenished by calling the source function
timeStepInterval = 1 # Timesteps between invocations of the source function

# The number of ions lost between source-function calls is:
ion_loss = ion_flux*timeStepInterval*ctrl.dt

# Compute a number-density increment to replace this loss
numberDensityIncrementHplus = ion_loss/(rmaxHplus-rminHplus)
# Check for positivity
if numberDensityIncrementHplus < 0:
    print "Check number density for species", ionSpeciesName, "is negative. Should be positive"
    sys.exit()

# Compute a value for thermalSpeed from the temperature in eV
temperature_eV = 0.1 # About 1000K

temp_joule = temperature_eV*MyPlasmaUnits_C.elem_charge
thermalSpeed = np_m.sqrt(2.0*temp_joule/ionMass)

driftVelocityHplus = (vdrift_1, vdrift_2, vdrift_3)

# Set desired particles-per-cell
numberPerCell = 1

## PS.3.1.3. Collect the parameters for this source in a dictionary
HplusSourceParams = {'species_name': ionSpeciesName,
#                     'source_distribution_type': sourceDistributionType,
                     'number_density': numberDensityIncrementHplus,
                     'thermal_speed': thermalSpeed,
                     'velocity_coordinate_system': velocity_coordinate_system,
                     'drift_velocity': driftVelocityHplus,
                     'timestep_interval': timeStepInterval,
                     'number_per_cell': numberPerCell}

if emitInput is True:
    print ""
    print "********** Hplus Source Parameters (HplusSourceParams) **********"
#    print "HplusSourceParams =", HplusSourceParams
    for k, v in HplusSourceParams.iteritems():
#        print "type(v) is ", type(v)
        if type(v) is np_m.float64:
            print " %s = %.3g" % (k, v)
        else:
            print '', k, '=', v
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

##### PS.3.2. Electron source near r=0


## PS.3.2.1. The species created in this source
speciesName = 'electron'

# Check that this species has been defined
if speciesName in particle_P.species_names:
    charge = particle_P.charge[speciesName]
    mass = particle_P.mass[speciesName]
else:
    print "The species", speciesName, "has not been defined. See Particle Species Attributes for defined species."
    sys.exit()

## PS.3.2.2. Source geometry for electrons: same as Hplus
rmin = rminHplus
rmax = rmaxHplus
electronSourceRegion = RectangularRegion_C(particle_P.pmesh_M, rmin, rmax)

if emitInput is True:
    print ""
    print "********** Electron Source Region Attributes **********"
    print electronSourceRegion
#    print "electronSourceRegion =", electronSourceRegion.__dict__
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)


## PS.3.2.3. Define the electron distribution function & parameters

# Name the particle-creation function:
particleGenerator = particle_P.create_maxwellian_particles
# !!!The function could provide it's own type
#sourceDistributionType = 'functional'

# Compute a number-density rate from the charge-density rate

numberDensityIncrement = numberDensityIncrementHplus
# Check for positivity
if numberDensityIncrement < 0:
    print "Check number density for species", speciesName, "is negative. Should be positive"
    sys.exit()

thermalSpeed = scr.vthermal_e

# Set a drift velocity
driftVelocity = driftVelocityHplus

# Set desired particles-per-cell
numberPerCell = 10

## PS.3.2.4. Collect the parameters for this source in a dictionary
electronSourceParams = {'species_name': speciesName,
                        #                     'source_distribution_type': sourceDistributionType,
                        'number_density': numberDensityIncrement,
                        'thermal_speed': thermalSpeed,
                        'velocity_coordinate_system': velocity_coordinate_system,
                        'drift_velocity': driftVelocity,
                        'timestep_interval': timeStepInterval,
                        'number_per_cell': numberPerCell}

if emitInput is True:
    print ""
    print "********** Electron Source Parameters (electronSourceParams) **********"
#    print "electronSourceParams =", electronSourceParams
    for k, v in electronSourceParams.iteritems():
        if type(v) is np_m.float64:
            print " %s = %.3g" % (k, v)
        else:
            print '', k, '=', v
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)


##### PS.3.3. Put all the particle source data into a dictionary

# The dictionary keys are mnemonic source names.
# The dictionary values are lists of tuples giving the
#   a. the physical source parameters.
#   b. particle creation function and
#   c. source region

## Note: the dictionary keys here are source-region names, NOT species names!
particleSourceDict = {'electron_source': (electronSourceParams, particleGenerator, electronSourceRegion),
                      'Hplus_source': (HplusSourceParams, particleGenerator, HplusSourceRegion),
                      }

particle_P.particle_source_dict = particleSourceDict

if emitInput is True:
    print ""
    print "********** Particle Generators (particleSourceDict) **********"
#    print "particleSourceDict =", particleSourceDict
    for k, v in particleSourceDict.iteritems():
        print '', str(k) + ':', v[1]
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)


# Check the particle attributes selected for output
particle_P.check_particle_output_parameters(ctrl)


############################## TRJ. Particle Trajectories ##############################

########## TRJ.1. Create input for a particle trajectory object

# Use an input object to collect initialization data for the trajectory object
trajin = TrajectoryInput_C()

trajin.maxpoints = None # Set to None to get every point
trajin.extra_points = 1 # Set to 1 to make sure one boundary-crossing can be accommodated
                        # (e.g., when a particle hits an absorbing boundary. Set to a
                        # larger value if there are multiple boundary reflections.

# Specify which particle variables to save.  This has the
# form of a numpy dtype specification.
# TODO. Use real 1D spherical coordinates(?).
trajin.explicit_dict = {'names': ['step', 't', 'x', 'ux', 'crossings', 'Ex'], 'formats': [int, np_m.float32, np_m.float32, np_m.float32, int, np_m.float32]}

########## TRJ.2. Create the trajectory object and attach it to the particle object.

# Note: no trajectory storage is created until particles marked with TRAJECTORY_FLAG
# are encountered.
traj_T = Trajectory_C(trajin, ctrl, particle_P.explicit_species, particle_P.implicit_species, particle_P.neutral_species)

## Attach this to the particle object.
particle_P.traj_T = traj_T


############################## SRC. Source terms for the electric potential ##############################

########## SRC.1. Charge-density obtained from kinetic particles

##### SRC.1.1. Create number-density arrays on the particle mesh for each kinetic species
number_density_element_type = 'Lagrange'
number_density_element_degree = 1
number_density_field_type = 'scalar'

for s in particle_P.species_names:
    particle_P.number_density_dict[s] = Field_C(mesh_M=particle_P.pmesh_M,
                                               element_type=number_density_element_type,
                                               element_degree=number_density_element_degree,
                                               field_type=number_density_field_type)

## SRC.1.2. Create a charge-density array on the field mesh
assembledCharge_F = Field_C(mesh_M=mesh_M,
                            element_type=number_density_element_type,
                            element_degree=number_density_element_degree,
                            field_type=number_density_field_type)

########## SRC.2 Charge-density from functions

# Attach this to the particle object?
#particle_P.charge_density = chargeDensity_F


############################## FS. Field Solver ##############################

########## FS.1. The scalar electric potential

# Use 'unset' to have a natural BC.

##### FS.1.1. Set boundary values of the potential
phiVals = {'rmin': 'unset',
           'rmax': 0.0
          }

# phiVals = {'rmin':  0.0,
#            'rmax': -scr.confining_potential,
#           }
           
if emitInput is True:
    print ""
    print "********** Potential Boundary Values (phiVals) **********"
#    print "phiVals =", phiVals
    for k, v in phiVals.iteritems():
        if type(v) is np_m.float64:
            print " %s = %.3g" % (k, v)
        else:
            print '', k, '=', v
#        print '', k, '=', v
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)


## Save these in a boundary-condition dictionary
phiBCs = dict( (bnd, [fieldBoundaryDict[bnd], phiVals[bnd]]) for bnd in fieldBoundaryDict.keys())

##### FS.1.2. Properties of the electric potential trial function
phiElementType = 'Lagrange'
phiElementDegree = 1
phiFieldType = 'scalar'

##### FS.1.3. Storage for the electric potential solution
phi_F = Field_C(mesh_M=mesh_M,
                element_type=phiElementType,
                element_degree=phiElementDegree,
                field_type=phiFieldType)

if emitInput is True:
    print ""
    print "********** Electric Potential Function Attributes **********"
#    print "phi_F =", phi_F.__dict__
    print phi_F
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

########## FS.2. The electric field (gradient of scalar potential)

##### FS.2.1. Properties of the electric field
if phiElementDegree == 1:
    # For linear elements, grad(phi) is discontinuous across
    # elements. To represent this field, we need Discontinuous Galerkin
    # elements.
    electricFieldElementType = 'DG'
else:
    electricFieldElementType = 'Lagrange'

##### FS.2.2. Storage for the electric field
negElectricField_F = Field_C(mesh_M=mesh_M,
                             element_type=electricFieldElementType,
                             element_degree=phiElementDegree-1,
                             field_type='vector')

if emitInput is True:
    print ""
    print "********** Electric Field Function Attributes **********"
#    print "negElectricField_F =", negElectricField_F.__dict__
    print negElectricField_F
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

########## FS.3. The potential field solver
    
##### FS.3.1. The potential field solver algorithm parameters
linearSolver = 'lu'
preconditioner = None

##### FS.3.2. Create the potential field solver object
fieldBoundaryMarker = mesh_M.field_boundary_marker
if emitInput is True:
    print ""
    print "********** Calling UserPoissonSolve_C constructor **********"
potentialsolve = UserPoissonSolve1DS_C(phi_F,
                                       linearSolver, preconditioner,
                                       fieldBoundaryMarker, phiBCs,
                                       neg_electric_field=negElectricField_F)

#print "Initial negElectricField_F:", negElectricField_F.function_values.get_local()

if emitInput is True:
    print ""
    print "********** Poisson Solver Attributes **********"
#    print "potentialsolve =", potentialsolve.__dict__
    print potentialsolve
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)


############################## EXT. External Fields ##############################

########## EXT.1. Create storage for an external electric field

# Use the same element type and degree as in the field solver.
if ctrl.apply_external_electric_field:
    externalElectricField_F = Field_C(mesh_M=mesh_M,
                                      element_type=electricFieldElementType,
                                      element_degree=phiElementDegree-1,
                                      field_type='vector')
    
# Put non-zero external electric field values into the electron source region
    externalElectricField_F.set_gaussian_values(scr.external_electric_field_amplitude, domain=electronSourceRegion)

    print "externalElectricField values:", externalElectricField_F.function_values.get_local()

#    import matplotlib.pyplot as mplot_m
#    df_m.plot(externalElectricField_F.function, title="External E")
#    mplot_m.show()
#    sys.exit()

############################## HST. Histories ##############################

########## HST.1. History input

##### HST.1.1. Particle histories
particleHistoryInput = ParticleHistoryInput_C()

particleHistoryInput.max_points = None # Set to None to get every point

# Globals: 'count', 'charge', 'KE', 'momentum'
particleHistoryInput.scalar_histories = ['count',]
#histin.particle_scalar_histories = ['count', 'charge',]

#fieldHistoryInput = FieldHistoryInput_C()
##### HST.1.2. Global field histories
#fieldHistoryInput.scalar_histories = ()
#histin.scalar_histories = ('E-field-energy', 'max-E-field',)

##### HST.1.3. Point-probe field histories
#fieldHistoryInput.probe_histories = ()
# histin.vector_histories = (('E', (0.0)), ('phi', (0.0)),
#                                 ('E', (2.0)), ('phi', (2.0)), 
#                                 ('E', (4.0)), ('phi', (4.0)),
#                                )

##### HST.1.4. Numerical parameter histories
#numericalHistoryInput = NumericalHistoryInput_C()

# ('average-Debye-length', 'maximum-Debye-length', 'Courant-number',
#  'average-streaming-number', 'maximum-streaming-number')
#numericalHistoryInput.histories = ()

########## HST.2. Create the history objects
#output.histories = History_C(histin, ctrl, particle_P.species_names)
particle_P.histories = ParticleHistory_C(particleHistoryInput, ctrl, particle_P.history_function_dict, particle_P.species_names)

if emitInputOnly is True:
    infoMsg = "\n(DnT INFO) %s exiting because emitInputOnly is True"  % (fileName)
    sys.exit(infoMsg)


############################## CHK. Check the simulation setup ##############################

if particleBoundaryDict is not None and particle_P.pmesh_M.particle_boundary_marker is None:
    print "A particleBoundaryDict was set, but particle_P.pmesh_M.particle_boundary_marker is None"
    sys.exit()
    

############################## INIT. Initialization ##############################

########## INIT.1. Do an initial field solve without particles

if plotInitialFields is True:
    print ""
    print "********** Compute and plot initial fields without charge-density from particles **********"
plotTitle = fileName + " initial field (no particles)"
potentialsolve.solve_for_phi(plot_flag=plotInitialFields, plot_title=plotTitle)

# Open the field output files and write the initial fields (no particles)
if ctrl.field_output_interval is not None:
    PotentialFieldOutputFile = df_m.File(ctrl.potential_field_output_file)
    phi_F.function.rename("phi", "phi_label")
    PotentialFieldOutputFile << phi_F.function

    # !!! VTK cannot write a vector field on a 1D grid.
    # ElectricFieldOutputFile = df_m.File(ctrl.electric_field_output_file)
    # negElectricField_F.function.rename("E", "E_label")
    # ElectricFieldOutputFile << negElectricField_F.function
    

########## INIT.2. Recompute initial fields including initial particle spacecharge
##########         and random electric field

# Open the particle output file and write the header
if ctrl.particle_output_interval is not None:
    particle_P.initialize_particle_output_file(ctrl)

##### INIT.2.1. Recalculate initial electric field
if particle_P.initial_particles_dict is not None:
    printFlags = {}
    for s in particle_P.species_names:
        printFlags[s] = False
    particle_P.initialize_particles(printFlags)
    # Write the initial particles to the output file
    if ctrl.particle_output_file is not None:
        particle_P.write_particles_to_file(ctrl)

    ## Record initial trajectory data
    # This will capture the FIRST point on the trajectories of marked particles
    if particle_P.initial_particles_dict is not None:
        if particle_P.traj_T is not None:
            particle_P.record_trajectory_data(ctrl.timeloop_count, ctrl.time, neg_E_field=potentialsolve.neg_electric_field)

    ##### INIT.2.1.1. Recompute the initial electrostatic field with particle spacecharge
    ## Accumulate number-density from kinetic particles
    ## and sum into the total charge-density.
    for s in particles_P.species_names:
        particles_P.accumulate_number_density(s)
        q = particles_P.charge[s]
        assembledCharge_F.multiply_add(particle_P.number_density_dict[s], multiplier=q)        

# This duplicates the above loop!
#    particle_P.accumulate_charge_density_from_particles(assembledCharge_F)

    if plotInitialFields is True:
        print ""
        print "********** Plot recomputed initial fields with charge-density from particles **********"
    plotTitleInitial = fileName + "Initial field (with particles)"
#    potentialsolve.solve_for_phi(plot_flag=plotInitialFields, plot_title=plotTitle)
    potentialsolve.solve_for_phi(assembled_charge=assembledCharge_F, assembled_charge_factor=scr.spacecharge_factor, plot_flag=plotInitialFields, plot_title=plotTitleInitial)
    

    ##### INIT.2.1.2. Add in a random electric field
    if ctrl.apply_external_electric_field is not None:
        # Compute new random fields
        externalElectricField_F.set_gaussian_values(scr.external_electric_field_amplitude, domain=electronSourceRegion)
        # Add these to the computed electric field
        potentialsolve.neg_electric_field.multiply_add(externalElectricField_F)
        
    # Write out the initial fields (with contribution from initial particles and external fields)
    if ctrl.field_output_interval is not None:
        phi_F.function.rename("phi", "phi_label")
        PotentialFieldOutputFile << phi_F.function

    # !!! VTK cannot write a vector field on a 1D grid.
    # ElectricFieldOutputFile = df_m.File(ctrl.electric_field_output_file)
    # negElectricField_F.function.rename("E", "E_label")
    # ElectricFieldOutputFile << negElectricField_F.function
    
else:
    print "(DnT INFO) %s\n\tThere are no initial particles in this simulation" % fileName

##### INIT.2.2. Record initial history data
if particle_P.histories is not None:
    particle_P.record_history_data(ctrl.timeloop_count, ctrl.time)
    

############################## RUN. Integrate forward in time ##############################

print ""
print "********** Integrate for %d timesteps from step %d, time %.3g, starting with %d particles **********" % (ctrl.n_timesteps, ctrl.timeloop_count, ctrl.time, particle_P.get_total_particle_count())
#print ""

for istep in xrange(ctrl.n_timesteps):

    ########## RUN.0. Increment the time counters
    ctrl.timeloop_count += 1
    ctrl.time += ctrl.dt

    ########## RUN.1. Advance the particles one timestep
    print ""
    print "***** Advance %d particles to step %d, time %.3g *****" % (particle_P.get_total_particle_count(), ctrl.timeloop_count, ctrl.time)
    particle_P.move_particles_in_electrostatic_field(ctrl, potentialsolve.neg_electric_field)

    ## Record trajectory data for all marked particles
    if particle_P.traj_T is not None:
        if (ctrl.timeloop_count % particle_P.traj_T.skip == 0):
            particle_P.record_trajectory_data(ctrl.timeloop_count, ctrl.time, neg_E_field=potentialsolve.neg_electric_field)

    ########## RUN.2. Add particles from source regions
    # Note: if a new particle gets marked for a trajectory, the first datum is
    # recorded inside the particle generator.
    if particle_P.particle_source_dict is not None:
        particle_P.add_more_particles(ctrl, neg_E_field=potentialsolve.neg_electric_field)

    ########## RUN.3. Write a snapshot of the particles
    if ctrl.particle_output_interval is not None:
        if ctrl.timeloop_count % ctrl.particle_output_interval == 0:
            particle_P.write_particles_to_file(ctrl)

    ########## RUN.4. Update the electric charge density
    for s in particle_P.species_names:
        particle_P.accumulate_number_density(s)
        q = particle_P.charge[s]
        assembledCharge_F.multiply_add(particle_P.number_density_dict[s], multiplier=q)

# this repeats the above loop:
#    particle_P.accumulate_charge_density_from_particles(assembledCharge_F)

    ########## RUN.5. Update the electric field

    ##### RUN.5.1. Set field plotting parameters
    if plotFieldsEveryTimestep is True:
        plotFields = True
        plotTitle = os.path.basename(__file__) + ": n=%d " % (ctrl.timeloop_count)
        print ""
        print "********** Plot field at timestep %d **********" % ctrl.timeloop_count
    elif plotFinalFields is True and ctrl.timeloop_count == ctrl.n_timesteps:
        plotFields = True
        plotTitle = os.path.basename(__file__) + ": n=%d " % (ctrl.timeloop_count)
        print ""
        print "********** Plot final field at timestep %d **********" % ctrl.timeloop_count
    else:
        plotFields = False
        plotTitle = None

    ##### RUN.5.2. Compute the new potential and E-field

# Apply assembled_factor before call by calling a "multiply" function?

    potentialsolve.solve_for_phi(assembled_charge=assembledCharge_F, assembled_charge_factor=scr.spacecharge_factor, plot_flag=plotFields, plot_title=plotTitle)

    ##### RUN.5.3. Add an external electric field
    if ctrl.apply_external_electric_field is not None:
        # Compute new random fields
        externalElectricField_F.set_gaussian_values(scr.external_electric_field_amplitude, domain=electronSourceRegion)
        # Add the external field to the self-electric field
        potentialsolve.neg_electric_field.multiply_add(externalElectricField_F, multiplier=None)

#        print "negElectricField_F:", negElectricField_F.function_values.get_local()
        print "neg_electric_field:", potentialsolve.neg_electric_field.function_values.get_local()
        
        
    ##### RUN.5.4. Write a snapshot of the fields
    if ctrl.field_output_interval is not None:
        if ctrl.timeloop_count % ctrl.field_output_interval == 0:
            phi_F.function.rename("phi", "phi_label")
            PotentialFieldOutputFile << phi_F.function
    
    ########## RUN.6. Record history data for this timestep
    if particle_P.histories is not None:
        if (ctrl.timeloop_count % particle_P.histories.skip == 0):
            particle_P.record_history_data(ctrl.timeloop_count, ctrl.time)

    if pauseAfterTimestep is True:
        pauseAfterTimestep=ctrl.get_keyboard_input(fileName)
        if pauseAfterTimestep is False:
            # Turn off plotting on each timestep
            plotFieldsEveryTimestep = False
    # End of the timestep loop (istep)

print ""
print "********** Exited the time loop at step %d, time %.3g, with %d particles **********" % (ctrl.timeloop_count, ctrl.time, particle_P.get_total_particle_count())
#print ""

#print fncName, "Exited the time loop", self.ctrl.timeloop_count, "to reach time", self.ctrl.time


############################## FIN. Finish the simulation ##############################

##### FIN.1. Write the last particle data and close the output file.
if ctrl.particle_output_file is not None:
    print ""
    print "********** Write final particle data to %s and close the file **********" % (ctrl.particle_output_file)
#    print ""
    particle_P.write_particles_to_file(ctrl, close_flag=True)

##### FIN.2. Record the LAST point on the particle trajectories.
if particle_P.traj_T is not None:
    print ""
    print "********** Record final particle trajectory data **********"
#    print ""
    particle_P.record_trajectory_data(ctrl.timeloop_count, ctrl.time, neg_E_field=potentialsolve.neg_electric_field)

##### FIN.3. Plot the particle trajectories on the particle mesh.
if os.environ.get('DISPLAY') is not None and plotTrajectoriesOnMesh is True:
    print ""
    print "********** Plot trajectories on particle mesh at end of simulation **********"
    plotTitle = os.path.basename(__file__) + ": "
    mesh = particle_P.pmesh_M.mesh
    holdPlot = True # Set to True to stop the plot from disappearing.
    particle_P.traj_T.plot_trajectories_on_mesh(mesh, plotTitle, hold_plot=holdPlot) # Plots trajectory spatial coordinates on top of the particle mesh

##### FIN.4. Plot the particle trajectories in phase-space.
if os.environ.get('DISPLAY') is not None and plotTrajectoriesInPhaseSpace is True:
    print ""
    print "********** Plot trajectories in phase-space at end of simulation **********"
    particle_P.traj_T.plot()

##### FIN.5. Plot the recorded history data.
if os.environ.get('DISPLAY') is not None and plotTimeHistories is True:
    print ""
    print "********** Plot recorded history data at end of simulation **********"
    if particle_P.histories is not None:
        particle_P.histories.plot()

##### FIN.6. Write final log message
print ""
print "********** Simulation ended normally at end of timestep %d, time %.3g **********" % (ctrl.timeloop_count, ctrl.time)
