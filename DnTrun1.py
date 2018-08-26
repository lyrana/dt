# DnT code in Python

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import os
import math
import importlib as im_M

import numpy as np_m
#import dolfin as df_m

import DnT_Module as DnT_M

import RecordedData_Module as Rec_M

#import IO_Module as IO_M
import UserUnits_Module as U_M

#print "This is DOLFIN Version", df_m.__version__

# The initialization, since this is an initial-value problem

# Initialize things that don't depend on input.

# Select the unit system to be used for input parameters.
Convert = U_M.MyPlasmaUnits_C

# Create the simulation object
system = DnT_M.DnTsystem_C()

# Write input with python.

#
# USER INPUT SECTION
#

# Create an instance of a DnTcontrol object
ctrl = DnT_M.DnTcontrol_C(use_mpi=True)

# Specify numerical precision
#prec = DnTparam.precision = np_m.float32
ctrl.precision = np_m.float64

# Timestep
ctrl.dt = 0.1

# Number of timesteps
ctrl.nsteps = 10 # doesn't usually need to be passed down to functions

#user_mesh_input = "UserMesh-QuarterCircle.py"
#DnTrun.user_mesh_fileobj = IO_M.open_file(user_mesh_input, 'r')

# may not need this?
# DnTrun.mesh_type = 'FE'
# if DnTrun.mesh_type not in DnTrun.mesh_type_options:
#     print 'Mesh type must be one of ', DnTrun.mesh_type_options
#     error_msg = "Mesh type must be one of %s" % DnTrun.mesh_type_options
#     sys.exit(error_msg)

# The particles live inside the spatial domain defined by the boundary
# of the mesh, so create the mesh before creating the particles
#

#
# GEOMETRY and MESH
#

# For this calculation, we use a finite-element mesh made with Dolfin

# The user can specify some geometry and mesh parameters here.

umi = DnT_M.DnTmeshInput_C()

# radial
umi.rmin, umi.rmax = 1.0, 5.0 # Mesh goes from rmin to rmax in radius
umi.nr = 10 # Number of divisions in r direction
umi.stretch = 1.3 # Stretch parameter

# theta, starts at 0
umi.tmax = math.pi/2 # theta extent of the mesh
umi.nt = 20  # Number of divisions in theta direction

# Specify the module where the user creates the geometry and mesh
# (i.e. the UserMesh_C class), and import it.
user_mesh_module = "UserMesh_y_Fields_FE2D_Module"
UMsh_M = im_M.import_module(user_mesh_module)

# Create the mesh
system.mesh_M = UMsh_M.UserMesh_C(umi)

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
system.charge_density = fI_M.Field_C(mesh=system.mesh_M,
                                       element_type=charge_density_element_type,
                                       element_degree=charge_density_element_degree,
                                       field_type=charge_density_field_type)

system.electrostatic_potential = fI_M.Field_C(mesh=system.mesh_M,
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

system.electric_field = fI_M.Field_C(mesh=system.mesh_M,
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
#psiCI = DnT_M.DnTpoissonSolveInput()
linear_solver = 'cg'
preconditioner = 'ilu'

#boundary_marker = 
phi_rmin = 0.0
phi_rmax = -1.5

computeEflag = True

system.poisson_solve = UFldSlv_M.UserPoissonSolve_C(linear_solver,
                                                    preconditioner,
                                                    phi_rmin,
                                                    phi_rmax,
                                                    computeEflag)


# DESCRIBE THE PARTICLES TO BE USED
ctrl.use_particles = True

if ctrl.use_particles == True:

    # Create an instance of the DnTparticleInput class
    pin = DnT_M.DnTparticleInput_C()

    pin.precision = ctrl.precision

    pin.force_components = ['x', 'y',]
    pin.force_precision = ctrl.precision

# Do these matter?  The mesh will dictate these? NO
# They matter because they're used to set the PARTICLE storage.


# Phase-space coordinates for particles
#    but the mesh spatial coords are given above.  Can't have particle coords that are different
#    from what the mesh holds? Yes, you can: Could have x, y on mesh and Fourier in z! Or a beam with x-y mesh and drifting in x,y,z external fields.

#    position_coordinates = ('x', 'y', 'z') # Usually, but not always, these are the same as the mesh coordinates, and have the same number of dimensions.

#if DnTrun.position_coordinates not in DnTrun.position_coordinate_options:
#    error_msg = "Position coordinates must be one of %s" % DnTrun.position_coordinate_options
#    sys.exit(error_msg)

# Could construct these names from the above: 'u' + comp

#    velocity_coordinates = ('ux', 'uy', 'uz')
#if DnTrun.velocity_coordinates not in DnTrun.velocity_coordinate_options:
#    error_msg = "Position coordinates must be one of %s" % DnTrun.velocity_coordinate_options
#    sys.exit(error_msg)

#    pin.phase_coords = position_coordinates + velocity_coordinates


# This could be moved into UserParticles, but it has to be the same for all species.
#    pin.position_coordinates = ('x', 'y', 'z')

# Define the particle species used
# need this?
#    DnTrun.user_particles_fileobj = IO_M.open_file(user_particles_input, 'r')

    # Give the properties of the particle species.  The charges and
    # masses are normally those of the physical plasma particles, and
    # not the computational macroparticles.  Macroparticle weights are
    # specified or computed in a separate file (see
    # user_particles_module below) giving the distribution functions,
    # and can vary from particle to particle.

    pin.particle_species = (('plasmaelectrons',
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
    pin.user_particles_module = "UserParticles_H_He_e"

    userParticlesModule = im_M.import_module(user_particles_module)
# Could use this variable if we need to import this module elsewhere:
    pin.user_particles_class = userParticlesClass = userParticlesModule.ParticleDistributions_C

#    ph = PH_M.Particles_Mesh(allSpecies)

# Length of the particle array segments
#    seg_len = DnTrun.PARTICLE_SEGMENT_LENGTH = 100

# Create the simulation object
#    run = DnT_M.DnTrun_C()


# Create the object that initializes and stores the particle species
# Don't need a seg_len here: Could be using storage that's not a segmented array.
# Pass DnTrun?, or storage info? Or just put seg_len in the SV class? Or a global?
#    DnTrun.particles = Part_M.Particles_C(particle_species, phase_coords, prec, seg_len, userParticlesClass, echoFlag=True)
    system.particle_P = Part_M.Particle_C(pin, printFlag=True)

# ??Generate the initial particle distributions
    rp_P = system.particle_P # abbrev for particle Class Instance
    printFlags = {}
    for sp in rp_P.names: printFlags[sp] = False
    # Turn plotting on for some species
    printFlags['plasmaelectrons'] = True

#    p_P.initialize_distributions(plotFlags)

# Interaction of fields and particles

    # Create the input for this...
    fpin = DnT_M.DnTfieldParticleInput_C()
    fpin.copy_field_mesh = False
    fpin.particle_integration_loop = 'loop-on-particles'
    fpin.force_components = ['x', 'y',]
    fpin.precision = ctrl.precision

    # ...and add the object to this run instance:
#    system.field_particle_interaction = FldPrt_M.Field_Particle_C(fpin, system.field_M, system.particle_P)

#    sys.exit()

# species type will be collected; e.g., allElectrons, allIons

# Create the initial state of the kinetic plasma from the particles and mesh
# DnTrun.particles_mesh = PrtMsh_M.Particles_UserMesh_C(DnTrun.particles, DnTrun.mesh)
# DnTrun.initial_state = PrtMsh_M.Particles_UserMesh_C(DnTrun.particles, DnTrun.mesh)

#    DnTrun.initial_state = PrtMsh_M.Particles_UserMesh_C(DnTrun.particles, DnTrun.mesh)

#ph.create_initial_particles(mesh)

# Plot the initial state
# DnTrun.plot(fields, particles)

# sys.exit()


# Particle trajectories

    # Create input for trajectories
    trajin = DnT_M.DnTtrajectoryInput_C()

    trajin.maxpoints = 1000 # Set to None to get every point

#    trajin.species_names = rp_P.species_names
# since these are not editable, should they go into an arg list?
#    trajin.explicit_species = 
#    trajin.implicit_species = 

    # Variables to save.  This has the form of a numpy dtype specification.
    trajin.explicitDict = {'names': ['x', 'ux', 'Ex'], 'formats': [np_m.float32, np_m.float32]}
    trajin.implicitDict = {'names': ['x', 'ux', 'phi'], 'formats': [np_m.float32, np_m.float32]}

#    run.particle_P.traj_T = Rec_M.Trajectory_C(trajin, ctrl)
    system.particle_P.traj_T = Rec_M.Trajectory_C(trajin, ctrl, rp_P.explicit_species, rp_P.implicit_species)

# Compute the initial fields

system.field_M.compute_electrostatic_potential(plotFlag=True)

# 
# Initialize particle species individually so they can be inspected
# Create particles from the user's input spec
# Create storage for number-densities computed from particles

if ctrl.use_particles == True:
    rp_P = system.particle_P
    rf_M = system.field_M
    for sp in rp_P.names:
        init_dist_type = rp_P.initial_distribution_type[sp]
        if init_dist_type == 'listed':
            rp_P.create_from_list(sp, printFlag[sp])
        elif rp_P.initial_distribution_type == 'functional':
            rp_P.create_from_functions(sp, printFlag[sp])
        elif rp_P.initial_distribution_type == 'particle_file':
            rp_P.create_from_file(sp, printFlag[sp])
        else:
            error_msg = "Unknown initial_distribution_type ", rp_P.initial_distribution_type, " in Main for species ", sp
            sys.exit(error_msg)
        # Create storage for particle number-densities, now that we
        # know how many species there are.  The species
        # number-densities are Field_C class variables.  The densities
        # are computed on the particle mesh.  The total charge-density
        # must be transferred to the field mesh for Poisson's
        # equation.
        rf_M.number_density[sp] = rf_M.fieldStorageCI.create_scalarField()

#
# Integrate the equations for the system forward in time
#
system.time_integrate_in_uniform_fields(ctrl)

system.time_integrate_in_electrostatic_field(ctrl)

system.time_integrate(ctrl)
    
