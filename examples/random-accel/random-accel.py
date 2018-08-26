#!/usr/bin/env python3

__version__ = 0.1
__author__ = 'Copyright (C) 2018 L. D. Hughes'
#__all__ = []

"""Apply a random electric field to an electron.
   Plot the energy of the electron vs. time.

   There's no self-electric field in this calculation.

"""
#### Commentary
#
#  The key parameters are:
#
#      1. The rms amplitude of the applied electric field
#         ('randomExternalElectricFieldAmplitude'). This determines the kick a particle
#         gets on each timestep.
#
#      2. The timestep ('ctrl.dt'). This determines for how long the field is applied
#         to the particles on each timestep.
#
#      3. The number of timesteps ('ctrl.n_timesteps'). This determines for how long
#         the particles spend in the random field.
#

import sys
import os
import numpy
import numpy as np_m
#import importlib as im_m
#import unittest

import dolfin as df_m

from DT_Module import DTcontrol_C
from Dolfin_Module import Mesh_C
from Dolfin_Module import Field_C
from Particle_Module import *
from RecordedData_Module import *

from UserMesh_y_Fields_FE_XYZ_Module import *

from UserUnits_Module import MyPlasmaUnits_C


### Physical constants

m_e = 1.0*MyPlasmaUnits_C.electron_mass
q_e = 1.0*MyPlasmaUnits_C.elem_charge

#-# User input #-#

numberOfParticles = 1000
numberOfTrajectories = 0
numberOfTrajectories = min(numberOfTrajectories, numberOfParticles)

ctrl = DTcontrol_C()

# Timestepping
ctrl.n_timesteps = 100 # 100000
ctrl.dt = 4.0e-7 #

# Initialize time counters
ctrl.timeloop_count = 0
ctrl.time = 0.0

# Set random seed
ctrl.random_seed = 1
np_m.random.seed(ctrl.random_seed)

# Electric field control
randomExternalElectricFieldAmplitude = 0.1
# Set the switches to empty dictionaries here if and only if they'll get set when the species are
# defined below. The default is that all defined electric fields will be applied to all
# species.
#ctrl.apply_solved_electric_field = {}
#ctrl.apply_random_external_electric_field = {}

# Initial conditions for the electron(s)

# Used in: p0 = (x0, x0, ux0, weight0, bitflags0, cell_index, unique_ID, crossings)
x0=0.0; ux0=0.0
weight0 = 2.0 # number of electrons per macroparticle
bitflags0 = 0b0 # bit flags variable
#cell_index = Mesh_C.NO_CELL
cell_index = 0
crossings = 0 # a counter

#-# Particles #-#

# Create an instance of the ParticleInput class
pin = ParticleInput_C()
# Initialize particles
pin.precision = numpy.float64
pin.particle_integration_loop = 'loop-on-particles'
pin.position_coordinates = ['x',] # Determines the particle-storage dimensions.
pin.force_components = ['x',] # Usually determined by the spatial mesh dimension.
pin.force_precision = numpy.float64

# Create a particle species object

speciesName = 'electron'
charge = -1.0*q_e
mass = 1.0*m_e
dynamics = 'explicit'
electron_S = ParticleSpecies_C(speciesName, charge, mass, dynamics)

# Set force-application switches for this species
# This is shown here as an example only.
if ctrl.apply_solved_electric_field is not None:
    ctrl.apply_solved_electric_field[speciesName] = False
if ctrl.apply_random_external_electric_field is not None:
    ctrl.apply_random_external_electric_field[speciesName] = True


# Add this species to particle input
pin.particle_species = (electron_S,)

# Create a Particle_C object
particle_P = Particle_C(pin, print_flag=True)

#-# Particle Trajectories #-#

# Create input for a particle trajectory object
trajin = TrajectoryInput_C()

trajin.maxpoints = None # Set to None to get every point
trajin.extra_points = 1 # Set to 1 to make sure one boundary-crossing can be accommodated
                        # (e.g., when a particle hits an absorbing boundary. Set to a
                        # larger value if there are multiple boundary reflections.

# Specify which particle variables to save.  This has the
# form of a numpy dtype specification.
#trajin.explicit_dict = {'names': ['step', 't', 'x', 'ux', 'crossings', 'Ex'], 'formats': [int, np_m.float32, np_m.float32, np_m.float32, int, np_m.float32]}

#HERE
trajin.explicit_dict = {'names': ['step', 't', 'ux', 'Ex_ext'], 'formats': [int, np_m.float32, np_m.float32, np_m.float32, int, np_m.float32]}

# Create the trajectory object and attach it to the particle object.

# Note: no trajectory storage is created until particles marked with TRAJECTORY_FLAG
# are encountered.
traj_T = Trajectory_C(trajin, ctrl, particle_P.explicit_species, particle_P.implicit_species, particle_P.neutral_species)

# Attach this to the particle object.
particle_P.traj_T = traj_T

if numberOfTrajectories != 0:
    skip = numberOfParticles/numberOfTrajectories
# Store copies of the particle.
for i in range(numberOfParticles):
    # Turn ON trajectory flag.
    if numberOfTrajectories != 0 and i % skip == 0:
        bitflags = bitflags0 | Particle_C.TRAJECTORY_FLAG
    else:
        bitflags = bitflags0
    unique_ID = Particle_C.UNIQUE_ID_COUNTER; Particle_C.UNIQUE_ID_COUNTER += 1
    # Trim the number of coordinates here if needed to match "position_coordinates"
    # variable in ParticleInput_C
    p0 = (x0, x0, ux0, weight0, bitflags, cell_index, unique_ID, crossings)
    p, pFullIndex = particle_P.pseg_arr[speciesName].put(p0)
    # Create a trajectory for the particle.
    if p['bitflags'] & particle_P.TRAJECTORY_FLAG != 0:
        dynamicsType = particle_P.dynamics[speciesName]
        particle_P.traj_T.create_trajectory(speciesName, pFullIndex, dynamicsType)

#-# Mesh input for a 1-cell 1D mesh #-#

# Specify the mesh parameters
umi1D = UserMeshInput_C()
umi1D.pmin = df_m.Point(0.0)
umi1D.pmax = df_m.Point(1.0)
umi1D.cells_on_side = (1,) # Need the comma to indicate a tuple

# Create the 1D particle mesh and add to the Particle_C object
plotTitle = os.path.basename(__file__) + ": " + ": X mesh"
pmesh1D = UserMesh_C(umi1D, compute_dictionaries=True, compute_tree=True, plot_flag=False, plot_title=plotTitle)

particle_P.pmesh_M = pmesh1D

# Set the initial cell index containing particle 
# particle_P.compute_mesh_cell_indices()


#-# Electric field properties #-#

electricFieldElementDegree = 0
electricFieldElementType = 'DG'

#-# The external electric field #-#

# The particle mover uses the negative electric field.
randomExternalElectricField_F = Field_C(mesh_M=pmesh1D,
                                     element_type=electricFieldElementType,
                                     element_degree=electricFieldElementDegree,
                                     field_type='vector')

# This sets electric-field values over the mesh or subdomain.
randomExternalElectricField_F.set_values(randomExternalElectricFieldAmplitude)

# This sets a value over a mesh or subdomain.
#externalElectricField_F.set_gaussian_values(externalElectricFieldAmplitude)

#-# Histories #-#

# Particle histories
particleHistoryInput = ParticleHistoryInput_C()

particleHistoryInput.max_points = None # Set to None to get every point

# Globals: 'count', 'charge', 'KE', 'momentum'
particleHistoryInput.scalar_histories = ['count', 'KE']
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

########## HST.2. Create the history objects
#output.histories = History_C(histin, ctrl, particle_P.species_names)

# Attach particle histories to the particle object
particle_P.histories = ParticleHistory_C(particleHistoryInput, ctrl, particle_P.history_function_dict, particle_P.species_names)

# Record initial trajectory datum
particle_P.record_trajectory_data(ctrl.timeloop_count, ctrl.time, external_E_field=randomExternalElectricField_F)

# Record initial particle history data
particle_P.record_history_data(ctrl.timeloop_count, ctrl.time)

#-# Integrate forward in time #-#

for istep in range(ctrl.n_timesteps):

    ctrl.timeloop_count += 1
    ctrl.time += ctrl.dt

#    print ""
#    print "***** Advance %d particles to step %d, time %.3g *****" % (particle_P.get_total_particle_count(), ctrl.timeloop_count, ctrl.time)
    particle_P.move_particles_in_electrostatic_field(ctrl, external_E_field=randomExternalElectricField_F, accel_only=True)

#    externalElectricField_F.set_gaussian_values(externalElectricFieldAmplitude)
    
    particle_P.record_trajectory_data(ctrl.timeloop_count, ctrl.time, external_E_field=randomExternalElectricField_F)

    particle_P.record_history_data(ctrl.timeloop_count, ctrl.time)    

print("")
print("********** Exited the time loop at step %d, time %.3g, with %d particles **********" % (ctrl.timeloop_count, ctrl.time, particle_P.get_total_particle_count()))

# Record the LAST point on the particle trajectory(ies).
particle_P.record_trajectory_data(ctrl.timeloop_count, ctrl.time, external_E_field=randomExternalElectricField_F)

# Plot the particle trajectories in phase-space.
if numberOfTrajectories != 0:
    particle_P.traj_T.plot(plot_vs_t_only=True)

if particle_P.histories is not None:
        particle_P.histories.plot()
