# DT code in Python

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import os
import math
import importlib as im_M

import numpy as np_M
#import dolfin as df_M

import DT_Module as DT_M
#import Field_Module as Fld_M
#import Particle_Module as Prt_M
#import Field_Particle_Module as FldPrt_M

import Trajectory_Module as Traj_M

#import Particles-Mesh_Module as PrtMsh_M

#import IO_Module as IO_M
import UserUnits_Module as U_M

#print "This is DOLFIN Version", df_M.DOLFIN_VERSION_STRING

# The initialization, since this is an initial-value problem

# Initialize things that don't depend on input.

# Select the unit system to be used for input parameters.
Convert = U_M.MyPlasmaUnits_C

# Create the simulation object
systemCI = DT_M.DTsystem_C()

# Write input with python.

#
# USER INPUT SECTION
#

# Create an instance of a DTcontrol object
ctrlCI = DT_M.DTcontrol_C(use_mpi=True)

# Specify numerical precision
#prec = DTparam.precision = np_M.float32
ctrlCI.precision = np_M.float64

# Timestep
ctrlCI.dt = 0.1

# Number of timesteps
ctrlCI.nsteps = 10 # doesn't usually need to be passed down to functions

#user_mesh_input = "UserMesh-QuarterCircle.py"
#DTrun.user_mesh_fileobj = IO_M.open_file(user_mesh_input, 'r')

# may not need this?
# DTrun.mesh_type = 'FE'
# if DTrun.mesh_type not in DTrun.mesh_type_options:
#     print 'Mesh type must be one of ', DTrun.mesh_type_options
#     error_msg = "Mesh type must be one of %s" % DTrun.mesh_type_options
#     sys.exit(error_msg)

# The particles live inside the spatial domain defined by the boundary
# of the mesh, so create the mesh before creating the particles
#

#
# GEOMETRY and MESH
#

# For this calculation, we use a finite-element mesh made with Dolfin

# The user can specify some geometry and mesh parameters here.

miCI = DT_M.DTmeshInput_C()

# radial
miCI.rmin, miCI.rmax = 1.0, 5.0 # Mesh goes from rmin to rmax in radius
miCI.nr = 10 # Number of divisions in r direction
miCI.stretch = 1.3 # Stretch parameter

# theta, starts at 0
miCI.tmax = math.pi/2 # theta extent of the mesh
miCI.nt = 20  # Number of divisions in theta direction

# Specify the module where the user creates the geometry and mesh
# (i.e. the UserMesh_C class), and import it.
user_mesh_module = "UserMesh_y_Fields_FE2D_Module"
UMsh_M = im_M.import_module(user_mesh_module)

# Create the mesh
systemCI.meshCI = UMsh_M.UserMesh_C(miCI)

#
# FIELDS needed for the calculation
#

# For this calculation, we use the Dolfin framework for finite-element fields.
field_infrastructure_module = 'Dolfin_Module'
fI_M = im_M.import_module(fieldInfrastructureModule)

# Allocate storage for the charge density, electrostatic potential, electric field
charge_density_element_type = 'Lagrange'
charge_density_element_degree = 1
charge_density_field_type = 'scalar'
systemCI.charge_density = fI_M.Field_C(mesh=systemCI.meshCI,
                                       element_type=charge_density_element_type,
                                       element_degree=charge_density_element_degree,
                                       field_type=charge_density_field_type)

systemCI.electrostatic_potential = fI_M.Field_C(mesh=systemCI.meshCI,
                                                element_type=charge_density_element_type,
                                                element_degree=charge_density_element_degree,
                                                field_type=charge_density_field_type)

if charge_density_element_degree == 1:
    # For linear elements, grad(phi) is discontinuous across
    # elements. To get these values, we need Discontinuous Galerkin
    # elements.
    electric_field_element_type = "DG"
else:
    electric_field_element_type = "Lagrange"

systemCI.electric_field = fI_M.Field_C(mesh=systemCI.meshCI,
                                       element_type=electric_field_element_type,
                                       element_degree=charge_density_element_degree-1,
                                       field_type='vector')

#
# FIELD-SOLVE SPECIFICATION
#

# Specify the module where the user has defined the field solves.  May
# be the same as the user mesh module above.
user_fieldsolve_module = "UserMesh_y_Fields_FE2D_Module"

UFldSlv_M = im_M.import_module(fiCI.user_fieldsolve_module)

# poissonSolveInput class contains specialized input for the poissinsolver
#psiCI = DT_M.DTpoissonSolveInput()
linear_solver = 'cg'
preconditioner = 'ilu'

boundary_marker = 
phi_rmin = 0.0
phi_rmax = -1.5

computeEflag = True

systemCI.poisson_solveCI = UFldSlv_M.UserPoissonSolve_C(linear_solver,
                                                    preconditioner,
                                                    phi_rmin,
                                                    phi_rmax,
                                                    computeEflag)


# DESCRIBE THE PARTICLES TO BE USED
ctrlCI.use_particles = True

if ctrlCI.use_particles == True:

    # Create an instance of the DTparticleInput class
    pinCI = DT_M.DTparticleInput_C()

    pinCI.precision = ctrlCI.precision

    pinCI.force_components = ['x', 'y',]
    pinCI.force_precision = ctrlCI.precision

# Do these matter?  The mesh will dictate these? NO
# They matter because they're used to set the PARTICLE storage.


# Phase-space coordinates for particles
#    but the mesh spatial coords are given above.  Can't have particle coords that are different
#    from what the mesh holds? Yes, you can: Could have x, y on mesh and Fourier in z! Or a beam with x-y mesh and drifting in x,y,z external fields.

#    position_coordinates = ('x', 'y', 'z') # Usually, but not always, these are the same as the mesh coordinates, and have the same number of dimensions.

#if DTrun.position_coordinates not in DTrun.position_coordinate_options:
#    error_msg = "Position coordinates must be one of %s" % DTrun.position_coordinate_options
#    sys.exit(error_msg)

# Could construct these names from the above: 'u' + comp

#    velocity_coordinates = ('ux', 'uy', 'uz')
#if DTrun.velocity_coordinates not in DTrun.velocity_coordinate_options:
#    error_msg = "Position coordinates must be one of %s" % DTrun.velocity_coordinate_options
#    sys.exit(error_msg)

#    pinCI.phase_coords = position_coordinates + velocity_coordinates


# This could be moved into UserParticles, but it has to be the same for all species.
#    pinCI.position_coordinates = ('x', 'y', 'z')

# Define the particle species used
# need this?
#    DTrun.user_particles_fileobj = IO_M.open_file(user_particles_input, 'r')

    # Give the properties of the particle species.  The charges and
    # masses are normally those of the physical plasma particles, and
    # not the computational macroparticles.  Macroparticle weights are
    # specified or computed in a separate file (see
    # user_particles_module below) giving the distribution functions,
    # and can vary from particle to particle.

    pinCI.particle_species = (('plasmaelectrons',
                       {'initial_distribution_type' : 'listed',
                        'charge' : -1.0*Convert.elem_charge,
                        'mass' : 1.0*Convert.electron_mass,
                        'dynamics' : 'implicit',
                        }
                         ),
                        ('Hplus', 
                         {'initial_distribution_type' : 'listed',
                          'charge' : 1.0*Convert.elem_charge,
                          'mass' : 1.0*Convert.AMU,
                          'dynamics' : 'explicit',
                          }
                         ),
                        ('He', 
                         {'initial_distribution_type' : None,
# A neutral has no charge. Remove this, or use as a flag?
                          'charge' : None,
                          'mass' : 4.0*Convert.AMU,
                          'dynamics' : 'explicit',
                          }
                         ),
                        )

# Specify the module where the user has given the particle distributions:
    pinCI.user_particles_module = "UserParticles_H_He_e"

    UPrt_M = im_M.import_module(user_particles_module)
# Could use this variable if we need to import this module elsewhere:
    pinCI.user_particles_class = UPrt_C = UPrt_M.ParticleDistributions_C

#    ph = PH_M.Particles_Mesh(allSpecies)

# Length of the particle array segments
#    seg_len = DTrun.PARTICLE_SEGMENT_LENGTH = 100

# Create the simulation object
#    runCI = DT_M.DTrun_C()


# Create the object that initializes and stores the particle species
# Don't need a seg_len here: Could be using storage that's not a segmented array.
# Pass DTrun?, or storage info? Or just put seg_len in the SV class? Or a global?
#    DTrun.particles = Part_M.Particles_C(particle_species, phase_coords, prec, seg_len, UPrt_C, echoFlag=True)
    systemCI.particleCI = Part_M.Particle_C(pinCI, printFlag=True)

# ??Generate the initial particle distributions
    rpCI = systemCI.particleCI # abbrev for particle Class Instance
    printFlags = {}
    for sp in rpCI.names: printFlags[sp] = False
    # Turn plotting on for some species
    printFlags['plasmaelectrons'] = True

#    pCI.initialize_distributions(plotFlags)

# Interaction of fields and particles

    # Create the input for this...
    fpinCI = DT_M.DTfieldParticleInput_C()
    fpinCI.copy_field_mesh = False
    fpinCI.particle_integration_loop = 'loop-on-particles'
    fpinCI.force_components = ['x', 'y',]
    fpinCI.precision = ctrlCI.precision

    # ...and add the object to this run instance:
#    systemCI.field_particle_interaction = FldPrt_M.Field_Particle_C(fpinCI, systemCI.fieldCI, systemCI.particleCI)

#    sys.exit()

# species type will be collected; e.g., allElectrons, allIons

# Create the initial state of the kinetic plasma from the particles and mesh
# DTrun.particles_mesh = PrtMsh_M.Particles_UserMesh_C(DTrun.particles, DTrun.mesh)
# DTrun.initial_state = PrtMsh_M.Particles_UserMesh_C(DTrun.particles, DTrun.mesh)

#    DTrun.initial_state = PrtMsh_M.Particles_UserMesh_C(DTrun.particles, DTrun.mesh)

#ph.create_initial_particles(mesh)

# Plot the initial state
# DTrun.plot(fields, particles)

# sys.exit()


# Particle trajectories

    # Create input for trajectories
    trajinCI = DT_M.DTtrajectoryInput_C()

    trajinCI.maxpoints = 1000 # Set to None to get every point

#    trajinCI.species_names = rpCI.species_names
# since these are not editable, should they go into an arg list?
#    trajinCI.explicit_species = 
#    trajinCI.implicit_species = 

    # Variables to save.  This has the form of a numpy dtype specification.
    trajinCI.explicitDict = {'names': ['x', 'ux', 'Ex'], 'formats': [np_M.float32, np_M.float32]}
    trajinCI.implicitDict = {'names': ['x', 'ux', 'phi'], 'formats': [np_M.float32, np_M.float32]}

#    runCI.particleCI.trajCI = Traj_M.Trajectory_C(trajinCI, ctrlCI)
    systemCI.particleCI.trajCI = Traj_M.Trajectory_C(trajinCI, ctrlCI, rpCI.explicit_species, rpCI.implicit_species)

# Compute the initial fields

systemCI.fieldCI.compute_electrostatic_potential(plotFlag=True)

# 
# Initialize particle species individually so they can be inspected
# Create particles from the user's input spec
# Create storage for number-densities computed from particles

if ctrlCI.use_particles == True:
    rpCI = systemCI.particleCI
    rfCI = systemCI.fieldCI
    for sp in rpCI.names:
        init_dist_type = rpCI.initial_distribution_type[sp]
        if init_dist_type == 'listed':
            rpCI.create_from_list(sp, printFlag[sp])
        elif rpCI.initial_distribution_type == 'functional':
            rpCI.create_from_functions(sp, printFlag[sp])
        elif rpCI.initial_distribution_type == 'particle_file':
            rpCI.create_from_file(sp, printFlag[sp])
        else:
            error_msg = "Unknown initial_distribution_type ", rpCI.initial_distribution_type, " in Main for species ", sp
            sys.exit(error_msg)
        # Create storage for particle number-densities, now that we
        # know how many species there are.  The species
        # number-densities are Field_C class variables.  The densities
        # are computed on the particle mesh.  The total charge-density
        # must be transferred to the field mesh for Poisson's
        # equation.
        rfCI.number_density[sp] = rfCI.fieldStorageCI.create_scalarField()

#
# Integrate the equations for the system forward in time
#
systemCI.time_integrate_in_uniform_fields(ctrlCI)

systemCI.time_integrate_in_electrostatic_field(ctrlCI)

systemCI.time_integrate(ctrlCI)
    
