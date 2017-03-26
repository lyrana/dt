# timestep, nsteps, 
# things you need for debugging, like procid

#__version__ = 0.1
#__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import numpy as np_M

class DTcontrol_C(object):
    """Top-level control of the run
    """

    def __init__(self, use_mpi = False):
        """ Initialize variables
        """

        # MPI procid?
        if use_mpi:
            self.use_mpi = True
            from mpi4py import MPI
            self.mpi_comm = MPI.COMM_WORLD
            self.mpi_rank = self.mpi_comm.Get_rank()
            self.mpi_size = self.mpi_comm.Get_size()

        self.precision = None

        # Timestep
        self.dt = None

        # Dimensions
#        self.particle_dimension = None

        self.use_particles = None
        self.PARTICLE_SEGMENT_LENGTH = None

        self.E0 = None
        self.B0 = None
        
        # Diagnostics

        return
#    def __init__(self, use_mpi = False):ENDDEF

#class DTcontrol_C(object):
    def write_variables():
         pass

#class DTcontrol_C(object):ENDCLASS

class DTsystem_C(object):
    """The physical system of interacting fields and particles.
    """

    def __init__(self):
        """ Initialize variables
        """

# May want things like this in order to call DT from a loop?
# or spawn off many runs?
# maybe don't need all of these:

        self.particleCI = None
        self.fieldCI = None

#        self.particles_mesh = None

        # the particle mesh is a copy of the field mesh
#        self.pmesh = df_M.Mesh(mesh)
        
#class DTsystem_C(object):
    def time_integrate_in_uniform_fields(self, ctrlCI):
        """Integrate the equations for the system forward in time, in
           fields that are uniform in space and time.
        """
        pCI = self.particleCI # abbrev for particle Class Instance

#        for istep in xrange(ctrlCI.nsteps):
#            for sp in pCI.names:
#                pCI.move_particles_in_uniform_fields(sp, ctrlCI)

        dt = ctrlCI.dt
        E0 = ctrlCI.E0
        for istep in xrange(ctrlCI.nsteps):
            # needs dt; doesn't need nsteps
            for sp in pCI.names:

                if self.pCI.get_species_particle_count(sp) == 0: continue # Skip if there are no particles in this species
                qmdt = self.qom[sp]*dt
                psa = pCI.pseg_arr[sp] # segmented array for this species

                psa.init_segment_loop()
                pseg = psa.get_next_segment()
                while isinstance(pseg, np_M.ndarray):
#                while pseg != None:
                    fpCI.compute_phi_at_particles(pseg)
                    pCI.move_particle_segment_in_uniform_fields(dt, qmdt, pseg, E0)
#                    print pseg['z']
                    pseg = psa.get_next_segment()

        return
#    def time_integrate_in_uniform_fields(self, ctrlCI):ENDDEF

#class DTsystem_C(object):
    def time_integrate_in_electrostatic_field(self, ctrlCI):
        """Integrate the equations for the system forward in time in
           an electrostatic field.
        """
        fCI = self.fields # abbrev for field Class Instance
        pCI = self.particleCI # abbrev for particle Class Instance
        fpCI = self.field_particle_interaction

#        Efield = fCI.fields.electric_field

        dt = ctrlCI.dt

        for istep in xrange(ctrlCI.nsteps):
            # needs dt; doesn't need nsteps

            # Gather particle trajectory data
            if pCI.trajCI.trajCI is not None:
                if istep % pCI.trajCI.skip == 0:
                    pCI.record_trajectory_data(fpCI)

            # Do the implicit species first
            if len(pCI.implicit_species != 0):
                self.iterate_implicit_electrostatic_particles(dt, pCI.implicit_species)

            # Then move the explicit species
            if len(pCI.explicit_species != 0):
                pCI.move_particles_in_electrostatic_field(dt, self.neg_electric_field)
#            self.step_particles_electrostatic_explicit(explicit_species)
                    
            # Apply boundary conditions on the particles
            
            fpCI.apply_particle_boundary_conditions()    
            
            # Compute the density
#            fpCI.compute_charge_density_from_particles()
            pCI.compute_charge_density(fCI)

            # Update the electrostatic potential
            fCI.compute_electrostatic_potential()

        return

            # Compute the electric potential

            # fCI.compute_electrostatic_potential()

            #     psa.init_segment_loop()
            #     pseg = psa.get_next_segment()


# Q: What's the dim of Eseg? same as the field (2D) or same as the particle coord (3D)
# A: same as the field


#position_coordinates = ('x', 'y', 'z') # Usually, but not always, these are the same as the mesh coordinates, and have the same number of dimensions.

# The names are 
#metadata = {'names' : pvars, 'formats': pvartypes}

#                Eseg = np.empty(len(pseg), dtype=item_dict))        


#            pCI.move_particles_in_electrostatic_field(dt, qmdt, Efield)
# P/C iterate on the electrons
# push electrons, resolve for phi



#                while pseg != None:
#                    E_at_particles = fpCI.compute_vector_field_at_particles(pseg, Efield)
#                    E_at_particles = fpCI.compute_vector_field_at_particles(pseg, Eseg)

#                    E_at_particles = fpCI.compute_E_at_particles(pseg, Eseg)
#                    E_at_particles = fpCI.compute_E_at_particles(pseg)


# why do the species & seg loop in here?  could push it to Field_Particles



#            pCI.move_implicit_species_in_electrostatic_field(dt, qmdt, fCI, fpCI)


# push explicit particles
#            pCI.move_explicit_species_in_electrostatic_field(dt, qmdt, fCI, fpCI)




#                    phi_at_particles = fpCI.compute_scalar_field_at_particles(pseg, phi)
#                    pCI.move_particle_segment_in_electrostatic_field(dt, qmdt, pseg, E)

#                    print pseg['z']
#                    pseg = psa.get_next_segment()

                # Maybe this guy is only for testing?
                # OR you build whatever calculation you want in this file!!!
#                pCI.move_particles_in_electrostatic_field(sp, ctrlCI)


 #       return
#    def time_integrate_in_electrostatic_field(self, ctrlCI):ENDDEF

#class DTsystem_C(object):
# just calls through to Particle_Module
    def step_particles_electrostatic_explicit(self, species_name):
        """Apply the standard explicit particle push in an
           electrostatic field.
        """
        pCI = self.particleCI
        fCI = self.fieldCI
        fpCI = self.fieldParticleCI

        # P/C loop on the field and particles
        
        for PCiter in range(PCiterations):
            # Push the implicit particles
#            fpCI.interpolate_potential_to_particles()
            pCI.move_particles_in_electrostatic_E(fpCI)
            # Compute the density
            fpCI.compute_density_from_particles()
            # Update the electrostatic potential
            fCI.compute_electrostatic_potential()

        return
#    def step_particles_electrostatic_explicit(self, species_name):ENDDEF


#class DTsystem_C(object):
    def iterate_implicit_electrostatic_particles(dt, implicit_species):
#    def step_particles_electrostatic_implicit(self, species_name):
        """The implicit push has a predictor/corrector
           (P/C) iteration of the particle push and the electrostatic
           potential.
        """
        pCI = self.particleCI
        fCI = self.fieldCI
        fpCI = self.fieldParticleCI

        # P/C loop on the field and particles
        
        for PCiter in range(PCiterations):
            # Push the implicit particles
#            fpCI.interpolate_potential_to_particles()
            pCI.move_particles_in_electrostatic_potential(dt, implicit_species, fCI, fpCI)
            # Compute the density of these species
            fpCI.compute_density_from_particles(implicit_species)
            # Update the electrostatic potential
            fCI.compute_electrostatic_potential()

        return
#    def iterate_implicit_electrostatic_particles(dt, implicit_species):ENDDEF

#class DTsystem_C(object):ENDCLASS

class DTparticleInput_C(object):
    """Particle input class.

       Contains the variables that describe the particles. The values are
       usually set by the user in MAIN.py.

       .. note:: The user can modify this class if different variables
                 are needed to specify the particles.

    """

    def __init__(self):

        # Usually set from ctrlCI.precision
        # Example: numpy.float64
        self.precision = None

        # Force components acting on the particles
        # e.g., ['x', 'y', 'z']
        self.force_components = None

        # Usually set from ctrlCI.precision
        # Example: numpy.float64
        self.force_precision = None

        # Values: 'loop-on-particles', 'loop-on-cells'
        self.particle_integration_loop = None

        # Determines the particle-storage dimensions
        # Example: ['x', 'y',]
        self.position_coordinates = None

# May want things like this in order to call DT from a loop?
# or spawn off many runs?
# maybe don't need all of these:
        self.particle_species = None

        # The initial particle mesh is a copy of the field mesh
#        self.pmesh = df_M.Mesh(mesh)
        self.pmeshCI = None

        # Module containing user-supplied particle distributions and
        # boundary-conditions.
        self.user_particles_module = None
        # The class containing distribution functions
        self.user_particles_class = None
        # The class containing particle boundary conditions
        self.user_particle_bcs_class = None

        return

#class DTparticleInput_C(object):ENDCLASS

class DTmeshInput_C(object):
    """Input for the field mesh
       The user can modify this for different mesh specifications.
       Field mesh, field solve control?  Could use to pass control things to the field mesh
    """

    def __init__(self):
        """ List the mesh variables that the user can set in MAIN.py
        """
#        self.mesh_type_options = ['FE', 'Cartesian']
#        self.mesh_type = None

        self.mesh_file = None

        self.user_mesh_input = None
        self.user_mesh_class = None

        self.precision = None
        self.mesh_class = None

        self.rmin = None
        self.rmax = None
        self.nr = None
        self.stretch = None
        
        self.tmax = None
        self.nt = None

        # Options: 'left, 'right', 'left/right', 'crossed'
        self.diagonal = None

        self.field_boundary_file = None
        # User-assigned names of mesh boundaries where Dirichlet
        # values are set.
        self.field_boundary_dict = None

        self.particle_boundary_file = None
        # User-assigned names of mesh boundaries where particle BCs
        # are set.
        self.particle_boundary_dict = None

        self.particle_source_file = None
        # User-assigned names of mesh regions where particles are
        # created
        self.particle_source_dict = None

# May want things like this in order to call DT from a loop?
# or spawn off many runs?
# maybe don't need all of these:
        self.meshCI = None
        self.pmeshCI = None

        # the particle mesh is a copy of the field mesh
#        self.pmeshCI = df_M.Mesh(meshCI)

        return
        
#class DTmeshInput_C(object):ENDCLASS

class DTpoissonSolveInput_C(object):
    """Input for the field solver(s).
       The user can modify this for different field solvers.
    """

    def __init__(self):
        """ List the field-solver parameters that the user
        can set in MAIN.py
        """

        self.user_poissonsolve_input = None
        self.user_poissonsolve_class = None

        self.meshCI = None

        self.element_type = None
        self.element_degree = None

        self.linear_solver = None
        self.preconditioner = None

        # Dirichlet BC object
        self.phi_BCs = None

        self.computeEflag = None

        return

#class DTpoissonSolveInput_C(object):ENDCLASS

class DTtrajectoryInput_C(object):
    """Parameters that specify particle trajectory
    """

    def __init__(self):
        """ Initialize variables
        """

        # Upper limit is the number of timesteps
        self.maxpoints = None

        self.npoints = None



        return

#class DTtrajectoryInput_C(object):ENDCLASS

class DToutput_C(object):
    """Output specifications
    """

    def __init__(self):
        """ Initialize variables
        """
        pass
    
# See Aleph input
        self.plot

        return

#class DToutput_C(object):ENDCLASS
