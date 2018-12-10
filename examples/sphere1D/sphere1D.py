#!/usr/bin/env python3

# 1D spherically-symmetric expansion of ions and electrons
# See sphere1D.ods for setting up parameters

__version__ = 0.1
__author__ = 'Copyright (C) 2017-2018 L. D. Hughes'
__all__ = ['Example_C',]

"""A spherical source region creates ion-electron pairs uniformly in space.  The
   electrons are hot and the ions are cold, but may be given a small-but-finite
   radial velocity.  The electrons are heated by an applied random electric
   field. The electrons pull the ions radially outward.

"""
#### Commentary

#  The radial extent of the simulation is from 'scr.rmin' to 'scr.rmax' (see SCR.1.)

#  The number of cells is 'scr.nr' (see SCR.1.)

#  The Debye length is specified by the number of cells in it (see SCR.1). This along
#  with the electron temperature yields the plasma density (see SCR.2).

#  The size of the source-region where ions and electrons are created is given as a number
#  of Debye lengths (see PS.3.1.2).

#  The ions are created with a radial drift-velocity 'vdrift_1' (see PS.3.1.3).

#  An ion flux out of the creation region is estimated from velocity 'vdrift_av'. This
#  is used to compute the number of ions to create per injection (see sphere1D.ods).


#  Key physical parameters are:

#       * The number of Debye-lengths in the simulation region.

#       * The electron temperature ('scr.T_e_eV'). This is used to compute the electron
#         density from the Debye length. The electron density and temperature determine
#         the electron pressure, and so the ambipolar electric field set up by the
#         electrons.

#       * 'ion_sound_speed_multiplier' (see SCR.3.). This sets the initial ion velocity, 'vdrift_1'.

#       * 'vdrift_av'. The user sets this to estimate the average ion drift speed in the
#         creation region, and this sets the rate at which new ions-electrons pairs are created
#         to replace the lost ions.

#       * The external random-electric-field amplitude
#         ('scr.random_external_electric_field_amplitude').  This is an external energy source
#         that heats the electrons in the source region.  Without this, cold electrons
#         will accumulate in this region while the hotter ones escape.

#       * The potential difference ('scr.confining_potential') between the center and the
#         edge of the simulation region.  This is the confining potential for the electrons.
#         (This is not used as a simulation parameter if the region extends to r = 0).

#  Key numerical parameters are:
#
#       * 'ctrl.n_timesteps' (CTRL.2.)
#       * 'numberPerCell' for the ions and electrons (PS.3.1.3.)
#         

#  Key output parameters are:
#
#       * 'ctrl.particle_output_interval'.
#

#### Section Index

### To search for a section: type a space after the last period of the section
### number below, e.g., PS.<space>, or PS.1.3.<space>

## IM. Imported Modules
## TO. Terminal Output flags
## PC. Physical Constants
## SCR. SCRatch-pad for calculating parameter values.
##     SCR.1. Calculate spatial parameters
##     SCR.2. Calculate electron parameters
##     SCR.3. Calculate ion parameters
##     SCR.4. Calculate field parameters
## CTRL. Simulation Control
##     CTRL.1. Titles
##     CTRL.2. Set timestep value and number of timesteps.
##     CTRL.3. Electric field control
##     CTRL.4. Output: Specify field and particle output files and output frequency.
## FM. Field Mesh
## PS. Particle Species
##     PS.1. Define the species, set the "apply-field" switches, and provide storage.
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
## DEN. Particle Densities
##      DEN.1. Create DoF number-density arrays on the particle mesh for each kinetic species
##      DEN.2. Create cell number-density arrays on the particle mesh for each kinetic species
## SRC. Source terms for the electric potential
##      SRC.1. Charge-density obtained from kinetic particles
##      SRC.2. Charge-density from functions
## FS. Field Solver
##     FS.1. The scalar electric potential
##           FS.1.1. Set boundary values of the potential
##     FS.2. The electric field (gradient of scalar potential)
##           FS.2.1. Properties of the electric field
##           FS.2.2. Storage for the electric field
##     FS.3. The potential field solver
##           FS.3.1. Compute dynamic electron permittivity (if implicit)
##           FS.3.2. Set the potential field solver algorithm parameters
##           FS.3.3. Create the potential field solver object
## EXT. External Fields
##      EXT.1. Create storage for an external electric field
## HST. Time-Histories
##      HST.1. History input
##             HST.1.1. Particle histories
##             HST.1.2. Global field histories
##             HST.1.3. Point-probe field histories
##             HST.1.4. Numerical parameter histories
##      HST.2. Create the history objects
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
##      RUN.4. Update the particle densities and electric charge density
##      RUN.5. Update the electric field
##             RUN.5.1. Set field-plotting parameters
##             RUN.5.2. Call the field solver
##             RUN.5.3. Add an external electric random field
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

import pseg_cpp
import dolfin_cpp

fileName = __file__+':'

############################## TO. Terminal output flags ##############################

pauseAfterWarning = False # If False, don't pause after warning messages.

emitInput = False # Write input values to standard output.
emitInputOnly = False # Write input values to standard output and then quit.
pauseAfterEmit = False # Wait for user input before continuing.  Can be changed to
                      # False by using ctrl.get_keyboard_input()

pauseAfterTimestep = False # Wait for user input before continuing to the next
                          # timestep.  Can be changed to False by using
                          # ctrl.get_keyboard_input()
                      
plotInitialFields = False
plotFieldsEveryTimestep = False # Plot phi and -E each timestep.
plotFinalFields = False

plotTrajectoriesOnMesh = False # Plot the trajectory particles on top of the mesh.
plotTrajectoriesInPhaseSpace = False # Plot the 2D phase-space of trajectory particles.

plotTimeHistories = False # Plot the recorded values vs. time.

############################## PC. Physical constants ##############################

# DO NOT MODIFY THESE VALUES!
# These are (abbreviations for) physical constants, and should not be modified.
# Values used in the simulation should be written in terms of these.
m_e = 1.0*MyPlasmaUnits_C.electron_mass
q_e = 1.0*MyPlasmaUnits_C.elem_charge

m_AMU = 1.0*MyPlasmaUnits_C.AMU

eps_0 = MyPlasmaUnits_C.epsilon_0
k_B = MyPlasmaUnits_C.boltzmann_constant


############################## SCR. Scratch-pad ##############################

# Use the scratch-pad to compute values needed to set up the simulation

if emitInput is True:
    print("")
    print("********** Scratch-pad for computing simulation parameters **********")

# Create an object to hold scratch-pad variables
scr = DTscratch_C()
# Create an object to hold simulation-control variables
ctrl = DTcontrol_C()

########## SCR.1. Calculate spatial parameters

# Start with the size of the computational domain
#scr.rmin = 1.0
scr.rmin = 0.0
scr.rmax = 8.0 # 20.0 40.0
# Set the spatial resolution
scr.nr = 50
if emitInput is True:
    print("The mesh extends from %.3g m to %.3g m in radius" % (scr.rmin, scr.rmax))
    print("There are %d cells in the mesh" % scr.nr)
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

# Compute a rough cell-size (The mesh can be stretched non-uniformly)
scr.dr_av = (scr.rmax-scr.rmin)/scr.nr
if emitInput is True:
    print("The average cell-size is %.3g m" % scr.dr_av)
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

# Specify the desired Debye length as a multiple of the cell-size
scr.lambda_D = 2.0*scr.dr_av
if emitInput is True:
    print("The electron Debye length is chosen to be %.3g m" % scr.lambda_D)
    print("There are roughly %.3g cells per electron Debye length" % (scr.lambda_D/scr.dr_av))
    print("There are roughly %.3g Debye lengths in the computational region" % ((scr.rmax-scr.rmin)/scr.lambda_D))
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

########## SCR.2. Calculate electron parameters

scr.electron_mass = m_e
# Specify the electron temperature
scr.T_e_eV = 0.02 # 6.0

# Set the electron streaming parameter, which sets the accuracy of the time-integration of
# electron orbits.
ctrl.electron_streaming_parameter = 0.3
    
# Convert this to an electron thermal velocity
scr.vthermal_e = np_m.sqrt(q_e*scr.T_e_eV/m_e)
if emitInput is True:
    print("T_e is set to %.3g eV, giving an electron thermal velocity of %.3g m/s" % (scr.T_e_eV, scr.vthermal_e))
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

# Compute the electron plasma frequency needed to get the above Debye length
scr.omega_e = scr.vthermal_e/scr.lambda_D
if emitInput is True:
    print("For the above electron Debye length, the electron plasma frequency is %.3g rad/s" % scr.omega_e)
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

# Compute the explicit timestep limit from the electron plasma frequency
scr.stable_dt_max = 2.0/scr.omega_e
if emitInput is True:
    print("This frequency gives a maximum stable timestep of %.3g s" % scr.stable_dt_max)
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

# Compute the timestep limit set by the desired streaming parameter
scr.streaming_dt_max = ctrl.electron_streaming_parameter*scr.dr_av/scr.vthermal_e
if emitInput is True:
    print("The electron streaming parameter %.3g puts a limit of %.3g s on the timestep" % (ctrl.electron_streaming_parameter, scr.streaming_dt_max))
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

# Compute the electron density from the plasma frequency
scr.electron_density = m_e*eps_0*(scr.omega_e/q_e)**2
if emitInput is True:
    print("The above electron plasma-frequency corresponds to an electron density of %.3g per m^{-3}" % scr.electron_density)
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

########## SCR.3. Calculate ion parameters

scr.ion_mass = 100.0*scr.electron_mass
if emitInput is True:
    print("The ratio of ion mass to electron mass is %.3g" % (scr.ion_mass/scr.electron_mass))
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)    

# Compute the ion sound speed
scr.ion_charge_state = 1.0
scr.ion_sound_speed = np_m.sqrt(scr.ion_charge_state*q_e*scr.T_e_eV/scr.ion_mass)
# If the mesh goes to r=0, the following is used to estimate the ion drift out of the source region:
scr.ion_sound_speed_multiplier = 0.1
if emitInput is True:
    print("The ion sound speed is %.3g m/s. Multiplier for initial ion drift is %.3g" % (scr.ion_sound_speed, scr.ion_sound_speed_multiplier))
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

#electron_plasma_frequency = q_e*np_m.sqrt(n_e/(m_e*eps_0))

########## SCR.4. Set and calculate field parameters

# Set the confining potential (in Volts)
#scr.confining_potential = 10.0*scr.T_e_eV

############################## CTRL. Simulation Control ##############################

##### CTRL.1. Titles

# Run identifier
ctrl.title = "sphere1D"
# Run author
ctrl.author = "tph"

##### CTRL.2. Timestepping

ctrl.n_timesteps = 100 # 1000

ctrl.dt = 1.0e-6 # sec
if ctrl.dt > scr.stable_dt_max:
    print("%s\n\tThe timestep %.3g exceeds stability limit by a factor of %.3g!" % (fileName, ctrl.dt, ctrl.dt/scr.stable_dt_max))
    if pauseAfterWarning is True:
        emitWarning=ctrl.get_keyboard_input(fileName)
elif ctrl.dt > scr.streaming_dt_max:
    print("%s\n\tThe timestep %.3g exceeds the streaming limit by a factor of %.3g!" % (fileName, ctrl.dt, ctrl.dt/scr.streaming_dt_max))
    if pauseAfterWarning is True:
        emitWarning=ctrl.get_keyboard_input(fileName)
#else:
#    print("%s\n\tTimestep is %.3g of stability limit" % (fileName, ctrl.dt/scr.stable_dt_max))
if emitInput is True:
    print("The timestep %.3g is %.3g times the stability limit and %.3g times the streaming limit" % (ctrl.dt, ctrl.dt/scr.stable_dt_max, ctrl.dt/scr.streaming_dt_max))
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

# Initialize time counters
ctrl.timeloop_count = 0
ctrl.time = 0.0

# Set random seed
ctrl.random_seed = 1
np_m.random.seed(ctrl.random_seed)

##### CTRL.3. Electric field control

# Electric field control

# Set the apply_field switches to empty dictionaries (i.e., {}) here if and only if you want to
# set per-species values in PS.1.1. below.  The default (if these dicts are set to None, or
# commented out) is that all defined electric fields will be applied to all
# species. If the following lines are commented out, then the forces WILL BE APPLIED
# to all species.

ctrl.apply_solved_electric_field = None
ctrl.apply_random_external_electric_field = None

spacechargeFactor = 1.0
randomExternalElectricFieldAmplitude = 0.0

##### CTRL.4. Output control

###Set parameters for writing field output files
ctrl.electron_density_output_file = ctrl.title + "-ne.xdmf"
ctrl.ion_density_output_file = ctrl.title + "-ni.xdmf"
ctrl.potential_field_output_file = ctrl.title + "-phi.xdmf"
# VTK output is deprecated: ctrl.potential_field_output_file = ctrl.title + "-phi.pvd"
# VTK output is deprecated: ctrl.electric_field_output_file = ctrl.title + "-E.pvd"

ctrl.cell_number_density_output_interval = 1
ctrl.field_output_interval = 1

### Set parameters for writing particle output to an h5part file
#ctrl.particle_output_file = "sphere1D.h5part"
ctrl.particle_output_file = ctrl.title + ".h5part"
ctrl.particle_output_interval = 1
# Available attributes: 'species_index' and any particle-record name
# e.g., 'weight', 'bit-flags', 'cell-index', 'unique_ID', 'crossings'
ctrl.particle_output_attributes = ('species_index', 'x', 'ux', 'weight', 'unique_ID', 'crossings')

if emitInput is True:
    print("")
    print("********** Simulation Control Attributes **********")
    print(ctrl)
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)


#    print ctrl.__dict__

#output = DToutput_C()


############################## FM. Field Mesh ##############################

########## FM.1. Provide parameters for a 1D mesh in the radial coordinate
#from UserMesh_y_Fields_Spherical1D_Module import *
umi = UserMeshInput1DS_C()

umi.rmin, umi.rmax = scr.rmin, scr.rmax # Mesh goes from rmin to rmax in radius
umi.nr = scr.nr # Number of divisions in r direction
umi.stretch = 1.0 # Stretch parameter

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
    print("")
    print("********** Field Boundary Dictionary (fieldBoundaryDict.__dict__) **********")
#    print "fieldBoundaryDict =", fieldBoundaryDict
    for k, v in fieldBoundaryDict.items():
        if type(v) is np_m.float64:
            print(" %s = %.3g" % (k, v))
        else:
            print('', k, '=', v)
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

############### PS.1. Define the species, set "apply-field" switches, and provide particle storage ###############

##### PS.1.1. Basic properties and storage
speciesName = 'electron'
charge = -1.0*q_e
mass = 1.0*m_e
dynamics = 'explicit'
electron_S = ParticleSpecies_C(speciesName, charge, mass, dynamics)

# Set force-application switches for this species
if ctrl.apply_solved_electric_field is not None:
    ctrl.apply_solved_electric_field[speciesName] = True
if ctrl.apply_random_external_electric_field is not None:
    ctrl.apply_random_external_electric_field[speciesName] = False

speciesName = 'Hplus'
charge = 1.0*q_e
mass = scr.ion_mass
dynamics = 'explicit'
Hplus_S = ParticleSpecies_C(speciesName, charge, mass, dynamics)

# Set force-application switches for this species
if ctrl.apply_solved_electric_field is not None:
    ctrl.apply_solved_electric_field[speciesName] = True
if ctrl.apply_random_external_electric_field is not None:
    ctrl.apply_random_external_electric_field[speciesName] = False

# Add the two species to particle input:
pin.particle_species = (electron_S, Hplus_S,)

if emitInput is True:
    print("")
    print("********** Particle Species Attributes (particle_species.__dict__) **********")
    for s in pin.particle_species:
        print("%s species:" % s.name)
        for k, v in s.__dict__.items():
            if type(v) is np_m.float64:
                print(" %s = %.3g" % (k, v))
            else:
                print('', k, '=', v)
#            print '', k, '=', v
        print("")
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

if emitInput is True:
    print("")
    print("********** Particle Input Attributes (pin.__dict__) **********")
#    print "Particle Input =", pin.__dict__
    for k, v in pin.__dict__.items():
            if type(v) is np_m.float64:
                print(" %s = %.3g" % (k, v))
            else:
                print('', k, '=', v)
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
    print("")
    print("********** Particle Boundary Dictionary (particleBoundaryDict) **********")
#    print "particleBoundaryDict =", particleBoundaryDict
    for k, v in particleBoundaryDict.items():
        if type(v) is np_m.float64:
            print(" %s = %.3g" % (k, v))
        else:
            print('', k, '=', v)
#        print '', k, '=', v
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

# Add the above dict to the mesh input
umi.particle_boundary_dict = particleBoundaryDict

if emitInput is True:
    print("")
    print("********** Field Mesh Input Attributes **********")
    print(umi)
#    print umi.__dict__
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

##### PS.2.2. Make the mesh and add it to the particle object

# First complete the mesh setup now that particle boundary-conditions have been
# set (see #FM1 above).
if emitInput is True:
    print("")
    print("********** Construct and plot the mesh **********")
plotTitle = os.path.basename(__file__) + ": "
mesh_M = UserMesh1DS_C(umi, compute_dictionaries=True, compute_tree=False, plot_flag=emitInput, plot_title=plotTitle)

# Write the mesh to an XML-format file so we can read it with a text-editor.
meshFile = df_m.File('sphere1D-mesh.xml')
meshFile << mesh_M.mesh

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
    errorMsg = "(DnT ERROR) %s The species %s has not been defined" % (fileName, ionSpeciesName)
    sys.exit(errorMsg)

## PS.3.1.2. Geometry of Hplus source region
rminHplus = scr.rmin
#rmaxHplus = rminHplus + 2.0*scr.dr_av
#rmaxHplus = rminHplus + 5.0*scr.lambda_D
rmaxHplus = 3.681 # Near exact boundary of cell if rmax=8
HplusSourceRegion = RectangularRegion_C(particle_P.pmesh_M, rminHplus, rmaxHplus)

if emitInput is True:
    print("")
    print("********** Hplus Source Region Attributes **********")
#    print "HplusSourceRegion =", HplusSourceRegion.__dict__
    print(HplusSourceRegion)
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

## PS.3.1.3. Define the Hplus generating function & parameters

# Name the particle-generating function:
particleGenerator = particle_P.create_maxwellian_particles
# !!!The function could provide it's own type
#sourceDistributionType = 'functional'

# For injection at some small inner radius, set the drift velocity >= the ion sound speed
velocity_coordinate_system = 'r_theta_z' # This sets the interpretation of the
                                         # following values
# vdrift_r = 2.64e4 # 1.1 times the sound speed
# vdrift_theta = 0.0
# vdrift_z = 0.0
#vdrift_1 = 1.1*scr.ion_sound_speed # 1.1 times the sound speed
#vdrift_1 = scr.ion_sound_speed_multiplier*scr.ion_sound_speed

# If the mesh goes to r=0, set the ion drift to zero ...
# vdrift_1 = 0.0
vdrift_2 = 0.0
vdrift_3 = 0.0
# ... and make an estimate of the average ion drift-velocity in the source region:
vdrift_av = scr.ion_sound_speed_multiplier*scr.ion_sound_speed

# 9jun18: Give the ions a speed to get them moving
vdrift_1 = vdrift_av

# Compute the flux of ions per unit area out of the source region (assuming the ion
# charge-density equals the electron charge-density):
ionFluxDensity = scr.electron_density*vdrift_av/scr.ion_charge_state

# The ions are replenished by calling the source function
timeStepInterval = 1 # Timesteps between invocations of the source function

# The number of ions lost per unit area between source-function calls is:
ionLossDensity = ionFluxDensity*timeStepInterval*ctrl.dt

# Compute a number-density increment to replace this loss
#sourceVolume = (4.0/3.0)*np_m.pi*(rmaxHplus**3 - rminHplus**3)
numberDensityIncrementHplus = 3.0*ionLossDensity/rmaxHplus
# Check for positivity
if numberDensityIncrementHplus < 0:
    print("Check number density for species", ionSpeciesName, "is negative. Should be positive")
    sys.exit()

# Compute a value for thermalSpeed from the temperature in eV
temperature_eV = 0.0 # 0.1 is about 1000K

temp_joule = temperature_eV*MyPlasmaUnits_C.elem_charge
thermalSpeed = np_m.sqrt(2.0*temp_joule/ionMass)

#driftVelocityHplus = (vdrift_1, vdrift_2, vdrift_3)
driftVelocityHplus = (vdrift_1,) # Truncate this velocity vector to match the
                                 # dimension of ParticleInput_C.position_coordinates.
                                 
# Extimate and write particle motion per timestep on the mesh
ionDriftCellsPerTimestep = driftVelocityHplus[0]*ctrl.dt/scr.dr_av
if emitInput is True:
    print("")
    print("Hplus ions drift %.3g average cells per timestep" % (ionDriftCellsPerTimestep))
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)
if ionDriftCellsPerTimestep >= 1.0:
    if pauseAfterWarning is True:
        emitWarning=ctrl.get_keyboard_input(fileName)


ionSpreadCellsPerTimestep = thermalSpeed*ctrl.dt/scr.dr_av
if emitInput is True:
    print("Hplus ions spread %.3g average cells per timestep" % (ionSpreadCellsPerTimestep))
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)    
if ionSpreadCellsPerTimestep >= 1.0:
    if pauseAfterWarning is True:
        emitWarning=ctrl.get_keyboard_input(fileName)

# Compute and output the ion density at edge of the creation region
numberDensityRate = numberDensityIncrementHplus/(timeStepInterval*ctrl.dt)
if emitInput is True:
    expectedDensity = (1.0/3.0)*numberDensityRate*rmaxHplus/driftVelocityHplus[0]
    print("Expected Hplus density at %.3g is %.3g" % (rmaxHplus, expectedDensity))
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)    
                                                   
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
    print("")
    print("********** Hplus Source Parameters (HplusSourceParams) **********")
#    print "HplusSourceParams =", HplusSourceParams
    for k, v in HplusSourceParams.items():
#        print "type(v) is ", type(v)
        if type(v) is np_m.float64:
            print(" %s = %.3g" % (k, v))
        else:
            print('', k, '=', v)
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

##### PS.3.2. Electron source near r=0


## PS.3.2.1. The species created in this source
speciesName = 'electron'

# Check that this species has been defined
if speciesName in particle_P.species_names:
    charge = particle_P.charge[speciesName]
    mass = particle_P.mass[speciesName]
else:
    print("The species", speciesName, "has not been defined. See Particle Species Attributes for defined species.")
    sys.exit()

## PS.3.2.2. Source geometry for electrons: same as Hplus
rmin = rminHplus
rmax = rmaxHplus
electronSourceRegion = RectangularRegion_C(particle_P.pmesh_M, rmin, rmax)

if emitInput is True:
    print("")
    print("********** Electron Source Region Attributes **********")
    print(electronSourceRegion)
#    print "electronSourceRegion =", electronSourceRegion.__dict__
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)


## PS.3.2.3. Define the electron distribution function & parameters

# Name the particle-creation function:
particleGenerator = particle_P.create_maxwellian_particles
# !!!The function could provide it's own type
#sourceDistributionType = 'functional'

# Add electrons at the same rate as Hplus ions (charge Z is same for both species in this case).

numberDensityIncrement = numberDensityIncrementHplus
# Check for positivity
if numberDensityIncrement < 0:
    print("Check number density for species", speciesName, "is negative. Should be positive")
    sys.exit()

thermalSpeed = scr.vthermal_e

# Set a drift velocity
driftVelocity = driftVelocityHplus

# Estimate and write out particle motion on the mesh per timestep
electronDriftCellsPerTimestep = driftVelocity[0]*ctrl.dt/scr.dr_av
if emitInput is True:
    print("")
    print("electrons drift %.3g average cells per timestep" % (electronDriftCellsPerTimestep))
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)
    
if electronDriftCellsPerTimestep >= 1.0:
    if pauseAfterWarning is True:
        emitWarning=ctrl.get_keyboard_input(fileName)

electronSpreadCellsPerTimestep = thermalSpeed*ctrl.dt/scr.dr_av
if emitInput is True:        
    print("electrons spread %.3g average cells per timestep" % (electronSpreadCellsPerTimestep))
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)
    
if electronSpreadCellsPerTimestep >= 1.0:
    if pauseAfterWarning is True:
        emitWarning=ctrl.get_keyboard_input(fileName)
    
# Set desired particles-per-cell
numberPerCell = 3 # 10

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
    print("")
    print("********** Electron Source Parameters (electronSourceParams) **********")
#    print "electronSourceParams =", electronSourceParams
    for k, v in electronSourceParams.items():
        if type(v) is np_m.float64:
            print(" %s = %.3g" % (k, v))
        else:
            print('', k, '=', v)
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
    print("")
    print("********** Particle Generators (particleSourceDict) **********")
#    print "particleSourceDict =", particleSourceDict
    for k, v in particleSourceDict.items():
        print('', str(k) + ':', v[1])
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

############################## DEN. Particle Densities ##############################

########## DEN.1. Create DoF number-density arrays on the particle mesh for each kinetic species

## Note: These arrays contain the inner product of the particle density function and
## the shape-function attached to each DoF.

dofNumberDensityElementType = 'Lagrange'
dofNumberDensityElementDegree = 1
dofNumberDensityFieldType = 'scalar'
dofNumberDensityDict_F = {s: None for s in particle_P.species_names}

for s in particle_P.species_names:
    dofNumberDensityDict_F[s] = Field_C(mesh_M=particle_P.pmesh_M,
                                        element_type=dofNumberDensityElementType,
                                        element_degree=dofNumberDensityElementDegree,
                                        field_type=dofNumberDensityFieldType)

########## DEN.2. Create cell number-density arrays on the particle mesh for each kinetic species

# These are used to compute the electron permittivity, and also as diagnostics.

cellNumberDensityElementType = 'DG'
cellNumberDensityElementDegree = 0
cellNumberDensityFieldType = 'scalar'
cellNumberDensityDict_F = {s: None for s in particle_P.species_names}

# Don't make the cell density an attribute of Particle_C (yet)
for s in particle_P.species_names:
    cellNumberDensityDict_F[s] = Field_C(mesh_M=particle_P.pmesh_M,
                                                     element_type=cellNumberDensityElementType,
                                                     element_degree=cellNumberDensityElementDegree,
                                                     field_type=cellNumberDensityFieldType)

############################## SRC. Source terms for the electric potential ##############################

########## SRC.1. Charge-density obtained from kinetic particles

# This is obained from the dof_number_density_dict[] arrays and so has same attributes
assembledCharge_F = Field_C(mesh_M=particle_P.pmesh_M,
                            element_type=dofNumberDensityElementType,
                            element_degree=dofNumberDensityElementDegree,
                            field_type=dofNumberDensityFieldType)

########## SRC.2. Charge-density from functions

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
    print("")
    print("********** Potential Boundary Values (phiVals) **********")
#    print "phiVals =", phiVals
    for k, v in phiVals.items():
        if type(v) is np_m.float64:
            print(" %s = %.3g" % (k, v))
        else:
            print('', k, '=', v)
#        print '', k, '=', v
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)


## Save these in a boundary-condition dictionary
phiBCs = dict( (bnd, [fieldBoundaryDict[bnd], phiVals[bnd]]) for bnd in list(fieldBoundaryDict.keys()))

##### FS.1.2. Properties of the electric potential trial function
phiElementType = 'Lagrange'
phiElementDegree = 1
phiFieldType = 'scalar'

##### FS.1.3. Storage for the electric potential solution
phi_F = Field_C(mesh_M=particle_P.pmesh_M,
                element_type=phiElementType,
                element_degree=phiElementDegree,
                field_type=phiFieldType)

if emitInput is True:
    print("")
    print("********** Electric Potential Function Attributes **********")
#    print "phi_F =", phi_F.__dict__
    print(phi_F)
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
electricFieldFieldType = 'vector'

##### FS.2.2. Storage for the electric field
negElectricField_F = Field_C(mesh_M=particle_P.pmesh_M,
                             element_type=electricFieldElementType,
                             element_degree=phiElementDegree-1,
                             field_type=electricFieldFieldType)

if emitInput is True:
    print("")
    print("********** Electric Field Function Attributes **********")
#    print "negElectricField_F =", negElectricField_F.__dict__
    print(negElectricField_F)
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

########## FS.3. The potential field solver
    
##### FS.3.1. Dynamic electron permittivity (if the electrons are treated implicitly)

# This is obtained from the cellNumberDensityDict_F[] arrays, and so has the same attributes.

electronPermittivity_F = Field_C(mesh_M=particle_P.pmesh_M,
                                 element_type=cellNumberDensityElementType,
                                 element_degree=cellNumberDensityElementDegree,
                                 field_type=cellNumberDensityFieldType)

##### FS.3.2. The potential field solver algorithm parameters
linearSolver = 'lu'
preconditioner = None

##### FS.3.3. Create the potential field solver object
fieldBoundaryMarker = particle_P.pmesh_M.field_boundary_marker
if emitInput is True:
    print("")
    print("********** Calling UserPoissonSolve_C constructor **********")
potentialsolve = UserPoissonSolve1DS_C(phi_F,
                                       linearSolver, preconditioner,
                                       fieldBoundaryMarker, phiBCs,
                                       neg_electric_field=negElectricField_F)

#print "Initial negElectricField_F:", negElectricField_F.function_values.get_local()

if emitInput is True:
    print("")
    print("********** Poisson Solver Attributes **********")
#    print "potentialsolve =", potentialsolve.__dict__
    print(potentialsolve)
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)


############################## EXT. External Fields ##############################

########## EXT.1. Create storage for an external electric field

# Use the same element type and degree as in the field solver.
if randomExternalElectricFieldAmplitude != 0.0:
    randomExternalElectricField_F = Field_C(mesh_M=particle_P.pmesh_M,
                                            element_type=electricFieldElementType,
                                            element_degree=phiElementDegree-1,
                                            field_type=electricFieldFieldType)
    # Put non-zero external electric field values into the electron source region
    randomExternalElectricField_F.set_values(randomExternalElectricFieldAmplitude, subdomain=electronSourceRegion)
    if emitInput is True:
        print("randomExternalElectricField values:", randomExternalElectricField_F.function_values.get_local())
        if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)        
else:
    randomExternalElectricField_F = None

#    import matplotlib.pyplot as mplot_m
#    df_m.plot(randomExternalElectricField_F.function, title="External E")
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
    print("A particleBoundaryDict was set, but particle_P.pmesh_M.particle_boundary_marker is None")
    sys.exit()

#TODO: add check on ctrl.apply_external_electric_field and apply_solved_electric_field, to ensure that the Field_C objects exist.    

############################## INIT. Initialization ##############################

########## INIT.1. Do an initial field solve without particles

if plotInitialFields is True:
    print("")
    print("********** Compute and plot initial fields without charge-density from particles **********")
plotTitle = fileName + " initial field (no particles)"
potentialsolve.solve_for_phi(plot_flag=plotInitialFields, plot_title=plotTitle)

# Open the field output files and write the initial fields (no particles)
if ctrl.field_output_interval is not None:
# VTK output is deprecated: potentialFieldOutputFile = df_m.File(ctrl.potential_field_output_file)
    potentialFieldOutputFile = ctrl.title + "-phi" + ".xdmf"
    potentialFieldOutputObj = df_m.XDMFFile(potentialFieldOutputFile)
    phi_F.function.rename("phi", "phi_label")
# VTK output is deprecated: potentialFieldOutputFile << phi_F.function
    potentialFieldOutputObj.write(phi_F.function, ctrl.time)

    # !!! VTK cannot write a vector field on a 1D grid.
    # ElectricFieldOutputFile = df_m.File(ctrl.electric_field_output_file)
    # negElectricField_F.function.rename("E", "E_label")
    # ElectricFieldOutputFile << negElectricField_F.function
    

########## INIT.2. Compute initial particle densities and recompute initial fields
##########         including particle spacecharge.

# Open XDMF output files for species cell number-densities
cellNumberDensityOutputObjs = {}
if ctrl.cell_number_density_output_interval is not None:
    for s in particle_P.species_names:
        cellNumberDensityOutputFile = ctrl.title + "-" + s + ".xdmf"
        cellNumberDensityOutputObjs[s] = df_m.XDMFFile(cellNumberDensityOutputFile)
        cellNumberDensityDict_F[s].function.rename(s+"-n", s + "_label")
    
# Open the particle output file and write the header
if ctrl.particle_output_interval is not None:
    particle_P.initialize_particle_output_file(ctrl)

##### INIT.2.1. Recalculate the initial electric field
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
            particle_P.record_trajectory_data(ctrl.timeloop_count, ctrl.time, neg_E_field=potentialsolve.neg_electric_field, external_E_field=randomExternalElectricField_F)

    ##### INIT.2.1.1. Recompute the initial electrostatic field with particle spacecharge

    ## Zero out the number-density and charge values before accumulation begins
    assembledCharge_F.set_values(0.0)
    for s in particle_P.species_names:
        dofNumberDensityDict_F[s].set_values(0.0)
        cellNumberDensityDict_F[s].set_values(0.0)
        
    ## Accumulate number-density from kinetic particles
    ## and sum into the total charge-density.
    for s in particle_P.species_names:
#        particle_P.accumulate_number_density(s, dofNumberDensityDict_F[s], cellNumberDensityDict_F[s])
        particle_P.accumulate_number_density_CPP(s, dofNumberDensityDict_F[s], cellNumberDensityDict_F[s])
        q = particle_P.charge[s]
        assembledCharge_F.multiply_add(dofNumberDensityDict_F[s], multiplier=q)        
# This duplicates the above loop!
#    particle_P.accumulate_charge_density_from_particles(assembledCharge_F)

    if plotInitialFields is True:
        print("")
        print("********** Plot recomputed initial fields with charge-density from particles **********")
    plotTitleInitial = fileName + "Initial field (with particles)"
#    potentialsolve.solve_for_phi(plot_flag=plotInitialFields, plot_title=plotTitle)
    potentialsolve.solve_for_phi(assembled_charge=assembledCharge_F, assembled_charge_factor=spacechargeFactor, plot_flag=plotInitialFields, plot_title=plotTitleInitial)
    
    ##### INIT.2.1.2. Add in an external electric field
#TODO: Change the following when code for a constant electric field is needed.
#      A random field is now applied inside the particle-advance loop.
#    if ctrl.apply_external_electric_field is not None:
#    if scr.external_electric_field_amplitude != 0.0:
        # Compute new random fields
#        randomExternalElectricField_F.set_gaussian_values(scr.random_external_electric_field_amplitude, domain=electronSourceRegion)
        # Add these to the computed (negative) electric field. 
#        potentialsolve.neg_electric_field.multiply_add(randomExternalElectricField_F, multiplier=-1.0)
        
    # Write out the initial fields (with contribution from initial particles and external fields)
    if ctrl.field_output_interval is not None:
# VTK output is deprecated: potentialFieldOutputFile << phi_F.function
        potentialFieldOutputObj.write(phi_F.function, ctrl.time)

    # !!! VTK cannot write a vector field on a 1D grid. VTK is no longer in Dolfin.
    # ElectricFieldOutputFile = df_m.File(ctrl.electric_field_output_file)
    # negElectricField_F.function.rename("E", "E_label")
    # ElectricFieldOutputFile << negElectricField_F.function
    
else:
    if emitInput is True:
        print("(DnT INFO) %s\n\tThere are no initial particles in this simulation" % fileName)
        if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

##### INIT.2.2. Record initial history data
if particle_P.histories is not None:
    particle_P.record_history_data(ctrl.timeloop_count, ctrl.time)
    

############################## RUN. Integrate forward in time ##############################

if emitInput is True:
    print("")
    print("********** Integrate for %d timesteps from step %d, time %.3g, starting with %d particles **********" % (ctrl.n_timesteps, ctrl.timeloop_count, ctrl.time, particle_P.get_total_particle_count()))
    if pauseAfterEmit is True: pauseAfterEmit=ctrl.get_keyboard_input(fileName)

for istep in range(ctrl.n_timesteps):

    ########## RUN.0. Increment the time counters
    ctrl.timeloop_count += 1
    ctrl.time += ctrl.dt

    ########## RUN.1. Advance the particles one timestep
    print("")
    print("***** Advance %d particles to step %d, time %.3g *****" % (particle_P.get_total_particle_count(), ctrl.timeloop_count, ctrl.time))
    particle_P.move_particles_in_electrostatic_field(ctrl, neg_E_field=potentialsolve.neg_electric_field, external_E_field=randomExternalElectricField_F)

    ## Record trajectory data for all marked particles
    if particle_P.traj_T is not None:
        if (ctrl.timeloop_count % particle_P.traj_T.skip == 0):
            particle_P.record_trajectory_data(ctrl.timeloop_count, ctrl.time, neg_E_field=potentialsolve.neg_electric_field, external_E_field=randomExternalElectricField_F)

    ########## RUN.2. Add particles from source regions
    # Note: if a new particle gets marked for a trajectory, the first datum is
    # recorded inside the particle generator.
    if particle_P.particle_source_dict is not None:
        particle_P.add_more_particles(ctrl, neg_E_field=potentialsolve.neg_electric_field)

    ########## RUN.3. Write a snapshot of the particles
    if ctrl.particle_output_interval is not None:
        if ctrl.timeloop_count % ctrl.particle_output_interval == 0:
            particle_P.write_particles_to_file(ctrl)

    ########## RUN.4. Update the particle densities and electric charge density
    ## Zero out the number-density and charge values before accumulation begins
    assembledCharge_F.set_values(0.0)
    for s in particle_P.species_names:
        dofNumberDensityDict_F[s].set_values(0.0)
        cellNumberDensityDict_F[s].set_values(0.0)

    for s in particle_P.species_names:
#        particle_P.accumulate_number_density(s, dofNumberDensityDict_F[s], cellNumberDensityDict_F[s])
        particle_P.accumulate_number_density_CPP(s, dofNumberDensityDict_F[s], cellNumberDensityDict_F[s])
        # Write a snapshot of the cell number densities
        if ctrl.cell_number_density_output_interval is not None:
            if ctrl.timeloop_count % ctrl.cell_number_density_output_interval == 0:
#                cellNumberDensity = dofNumberDensity_F.function.vector().get_local()
#                print("species", s, "has cell density", cellNumberDensity
                cellNumberDensityOutputObjs[s].write(cellNumberDensityDict_F[s].function, ctrl.time)
        q = particle_P.charge[s]
        assembledCharge_F.multiply_add(dofNumberDensityDict_F[s], multiplier=q)

# this repeats the above loop:
#    particle_P.accumulate_charge_density_from_particles(assembledCharge_F)

    ########## RUN.5. Update the electric field

    ##### RUN.5.1. Set field plotting parameters
    if plotFieldsEveryTimestep is True:
        plotFields = True
        plotTitle = os.path.basename(__file__) + ": n=%d " % (ctrl.timeloop_count)
        print("")
        print("********** Plot field at timestep %d **********" % ctrl.timeloop_count)
    elif plotFinalFields is True and ctrl.timeloop_count == ctrl.n_timesteps:
        plotFields = True
        plotTitle = os.path.basename(__file__) + ": n=%d " % (ctrl.timeloop_count)
        print("")
        print("********** Plot final field at timestep %d **********" % ctrl.timeloop_count)
    else:
        plotFields = False
        plotTitle = None

    ##### RUN.5.2. Compute the new potential and E-field

# Apply assembled_factor before call by calling a "multiply" function?

    potentialsolve.solve_for_phi(assembled_charge=assembledCharge_F, assembled_charge_factor=spacechargeFactor, plot_flag=plotFields, plot_title=plotTitle)

    ##### RUN.5.3. Add an external electric field
# TODO: Modify the following for a non-random external field, when needed.    
#    if ctrl.apply_external_electric_field is not None:
        # Compute new random fields
#        randomExternalElectricField_F.set_gaussian_values(scr.random_external_electric_field_amplitude, domain=electronSourceRegion)
        # Add the external field to the self-electric field
#        potentialsolve.neg_electric_field.multiply_add(randomExternalElectricField_F, multiplier=None)

#        print "neg_electric_field:", potentialsolve.neg_electric_field.function_values.get_local()


    ##### RUN.5.4. Write a snapshot of the fields
    if ctrl.field_output_interval is not None:
        if ctrl.timeloop_count % ctrl.field_output_interval == 0:
# VTK output is deprecated: potentialFieldOutputFile << phi_F.function
            potentialFieldOutputObj.write(phi_F.function, ctrl.time)            
    
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

print("")
print("********** Exited the time loop at step %d, time %.3g, with %d particles **********" % (ctrl.timeloop_count, ctrl.time, particle_P.get_total_particle_count()))
#print ""

#print fncName, "Exited the time loop", self.ctrl.timeloop_count, "to reach time", self.ctrl.time


############################## FIN. Finish the simulation ##############################

##### FIN.1. Write the last particle data and close the output file.
if ctrl.particle_output_file is not None:
    print("")
    print("********** Write final particle data to %s and close the file **********" % (ctrl.particle_output_file))
#    print ""
    particle_P.write_particles_to_file(ctrl, close_flag=True)

##### FIN.2. Record the LAST point on the particle trajectories.
if particle_P.traj_T is not None:
    print("")
    print("********** Record final particle trajectory data **********")
#    print ""
    particle_P.record_trajectory_data(ctrl.timeloop_count, ctrl.time, neg_E_field=potentialsolve.neg_electric_field, external_E_field=randomExternalElectricField_F)

##### FIN.3. Plot the particle trajectories on the particle mesh.
if os.environ.get('DISPLAY') is not None and plotTrajectoriesOnMesh is True:
    print("")
    print("********** Plot trajectories on particle mesh at end of simulation **********")
    plotTitle = os.path.basename(__file__) + ": "
    mesh = particle_P.pmesh_M.mesh
    holdPlot = True # Set to True to stop the plot from disappearing.
    particle_P.traj_T.plot_trajectories_on_mesh(mesh, plotTitle, hold_plot=holdPlot) # Plots trajectory spatial coordinates on top of the particle mesh

##### FIN.4. Plot the particle trajectories in phase-space.
if os.environ.get('DISPLAY') is not None and plotTrajectoriesInPhaseSpace is True:
    print("")
    print("********** Plot trajectories in phase-space at end of simulation **********")
    particle_P.traj_T.plot()

##### FIN.5. Plot the recorded history data.
if os.environ.get('DISPLAY') is not None and plotTimeHistories is True:
    print("")
    print("********** Plot recorded history data at end of simulation **********")
    if particle_P.histories is not None:
        particle_P.histories.plot()

##### FIN.6. Write final log message
print("")
print("********** Simulation ended normally at end of timestep %d, time %.3g **********" % (ctrl.timeloop_count, ctrl.time))
