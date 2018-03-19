# timestep, n_timesteps, 
# things you need for debugging, like procid

#__version__ = 0.1
#__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import numpy as np_m

#STARTCLASS
class DTcontrol_C(object):
    """A class containing user input to control the numerical simulation.
    """

    def __init__(self, use_mpi = False):
        """ Initialize variables
        """

        # Known coordinate systems
        self.coordinate_systems = ('Cartesian', '1D-spherical-radius')

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

        # Spatial coordinates
        self.coordinate_system = None

        # Dimensions
#        self.particle_dimension = None

        self.use_particles = None
        self.PARTICLE_SEGMENT_LENGTH = None

        self.E0 = None
        self.B0 = None
        
        # Random number seed
        self.random_seed = np_m.random.seed(1)

        # Diagnostics

        self.potential_field_output_file = None
        self.electric_field_output_file = None
        self.field_output_interval = None

        self.particle_output_file = None
        self.particle_output_interval = None
        self.particle_output_attributes = None
#        self.particle_output_attribute_dtypes = None

        # Run information

        self.timeloop_count = None
        self.time = None

        return
#    def __init__(self, use_mpi = False):ENDDEF

#class DTcontrol_C(object):
    def __str__(self):
        """Print the class members.
        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():'
        print fncName

        print_string = ' title = ' + str(self.title)
        print_string += '\n author = ' + str(self.author)
        print_string += '\n dt = ' + str(self.dt)
        print_string += '\n n_timesteps = ' + str(self.n_timesteps)
        print_string += '\n particle_output_file = ' + str(self.particle_output_file)
        print_string += '\n particle_output_attributes = ' + str(self.particle_output_attributes)

        return print_string
#     def __str__(self):ENDDEF

#class DTcontrol_C(object):
    def get_keyboard_input(self, callerName):
        """
           Reads input from the keyboard.

           :param callerName: Name of the function calling this function
           :param pauseAfterEmit: Boolean to pause after emitting output

        """

        inChar = raw_input("Press Enter to continue, q to quit, ! to stop pausing...")
        if inChar == 'q':
            infoMsg = "\n(DnT INFO) %s Exiting because %s was entered"  % (callerName, inChar)
            sys.exit(infoMsg)
        elif inChar == '!':
            pauseAfterEmit = False
        else:
            pauseAfterEmit = True

        return pauseAfterEmit
#    def process_keyboard_input(self):ENDDEF

#class DTcontrol_C(object):ENDCLASS

#STARTCLASS
class DTsystem_C(object):
    """The physical system of interacting fields and particles.
    """

    def __init__(self):
        """ Initialize variables
        """

# May want things like this in order to call DT from a loop?
# or spawn off many runs?
# maybe don't need all of these:

        self.particle_P = None
        self.field_M = None

#        self.particles_mesh = None

        # the particle mesh is a copy of the field mesh
#        self.pmesh = df_m.Mesh(mesh)
        
#class DTsystem_C(object):
    def time_integrate_in_uniform_fields(self, ctrl):
        """Integrate the equations for the system forward in time, in
           fields that are uniform in space and time.
        """
        p_P = self.particle_P # abbrev for particle Class Instance

#        for istep in xrange(ctrl.n_timesteps):
#            for sp in p_P.names:
#                p_P.move_particles_in_uniform_fields(sp, ctrl)

        dt = ctrl.dt
        E0 = ctrl.E0
        for istep in xrange(ctrl.n_timesteps):
            # needs dt; doesn't need n_timesteps
            for sp in p_P.names:

                if self.p_P.get_species_particle_count(sp) == 0: continue # Skip if there are no particles in this species
                qmdt = self.qom[sp]*dt
                psa = p_P.pseg_arr[sp] # segmented array for this species

                psa.init_segment_loop()
                pseg = psa.get_next_segment()
                while isinstance(pseg, np_m.ndarray):
#                while pseg != None:
                    fpCI.compute_phi_at_particles(pseg)
                    p_P.move_particle_segment_in_uniform_fields(dt, qmdt, pseg, E0)
#                    print pseg['z']
                    pseg = psa.get_next_segment()

        return
#    def time_integrate_in_uniform_fields(self, ctrl):ENDDEF

#class DTsystem_C(object):
    def time_integrate_in_electrostatic_field(self, ctrl):
        """Integrate the equations for the system forward in time in
           an electrostatic field.
        """
        f_M = self.fields # abbrev for field Class Instance
        p_P = self.particle_P # abbrev for particle Class Instance
        fpCI = self.field_particle_interaction

#        Efield = f_M.fields.electric_field

        dt = ctrl.dt

        for istep in xrange(ctrl.n_timesteps):
            # needs dt; doesn't need n_timesteps

            # Gather particle trajectory data
            if p_P.traj_T.traj_T is not None:
                if istep % p_P.traj_T.skip == 0:
                    p_P.record_trajectory_data(fpCI)

            # Do the implicit species first
            if len(p_P.implicit_species != 0):
                self.iterate_implicit_electrostatic_particles(dt, p_P.implicit_species)

            # Then move the explicit species
            if len(p_P.explicit_species != 0):
                p_P.move_particles_in_electrostatic_field(dt, self.neg_electric_field)
#            self.step_particles_electrostatic_explicit(explicit_species)
                    
            # Apply boundary conditions on the particles
            
            fpCI.apply_particle_boundary_conditions()    
            
            # Compute the density
#            fpCI.compute_charge_density_from_particles()
            p_P.compute_charge_density(f_M)

            # Update the electrostatic potential
            f_M.compute_electrostatic_potential()

        return

            # Compute the electric potential

            # f_M.compute_electrostatic_potential()

            #     psa.init_segment_loop()
            #     pseg = psa.get_next_segment()


# Q: What's the dim of Eseg? same as the field (2D) or same as the particle coord (3D)
# A: same as the field


#position_coordinates = ('x', 'y', 'z') # Usually, but not always, these are the same as the mesh coordinates, and have the same number of dimensions.

# The names are 
#metadata = {'names' : pvars, 'formats': pvartypes}

#                Eseg = np.empty(len(pseg), dtype=item_dict))        


#            p_P.move_particles_in_electrostatic_field(dt, qmdt, Efield)
# P/C iterate on the electrons
# push electrons, resolve for phi



#                while pseg != None:
#                    E_at_particles = fpCI.compute_vector_field_at_particles(pseg, Efield)
#                    E_at_particles = fpCI.compute_vector_field_at_particles(pseg, Eseg)

#                    E_at_particles = fpCI.compute_E_at_particles(pseg, Eseg)
#                    E_at_particles = fpCI.compute_E_at_particles(pseg)


# why do the species & seg loop in here?  could push it to Field_Particles



#            p_P.move_implicit_species_in_electrostatic_field(dt, qmdt, f_M, fpCI)


# push explicit particles
#            p_P.move_explicit_species_in_electrostatic_field(dt, qmdt, f_M, fpCI)




#                    phi_at_particles = fpCI.compute_scalar_field_at_particles(pseg, phi)
#                    p_P.move_particle_segment_in_electrostatic_field(dt, qmdt, pseg, E)

#                    print pseg['z']
#                    pseg = psa.get_next_segment()

                # Maybe this guy is only for testing?
                # OR you build whatever calculation you want in this file!!!
#                p_P.move_particles_in_electrostatic_field(sp, ctrl)


 #       return
#    def time_integrate_in_electrostatic_field(self, ctrl):ENDDEF

#class DTsystem_C(object):
# just calls through to Particle_Module
    def step_particles_electrostatic_explicit(self, species_name):
        """Apply the standard explicit particle push in an
           electrostatic field.
        """
        p_P = self.particle_P
        f_M = self.field_M
        fpCI = self.fieldParticleCI

        # P/C loop on the field and particles
        
        for PCiter in range(PCiterations):
            # Push the implicit particles
#            fpCI.interpolate_potential_to_particles()
            p_P.move_particles_in_electrostatic_E(fpCI)
            # Compute the density
            fpCI.compute_density_from_particles()
            # Update the electrostatic potential
            f_M.compute_electrostatic_potential()

        return
#    def step_particles_electrostatic_explicit(self, species_name):ENDDEF


#class DTsystem_C(object):
    def iterate_implicit_electrostatic_particles(dt, implicit_species):
#    def step_particles_electrostatic_implicit(self, species_name):
        """The implicit push has a predictor/corrector
           (P/C) iteration of the particle push and the electrostatic
           potential.
        """
        p_P = self.particle_P
        f_M = self.field_M
        fpCI = self.fieldParticleCI

        # P/C loop on the field and particles
        
        for PCiter in range(PCiterations):
            # Push the implicit particles
#            fpCI.interpolate_potential_to_particles()
            p_P.move_particles_in_electrostatic_potential(dt, implicit_species, f_M, fpCI)
            # Compute the density of these species
            fpCI.compute_density_from_particles(implicit_species)
            # Update the electrostatic potential
            f_M.compute_electrostatic_potential()

        return
#    def iterate_implicit_electrostatic_particles(dt, implicit_species):ENDDEF

#class DTsystem_C(object):ENDCLASS

#STARTCLASS
class DToutput_C(object):
    """Output specifications
    """

    def __init__(self):
        """ Initialize variables
        """
        self.histories = None
    
        return

#class DToutput_C(object):ENDCLASS

    
#STARTCLASS
class DTscratch_C(object):
    """Contains scratch-pad values
    """

    def __init__(self):
        """ Initialize variables
        """
    
        return

#class DTscratch_C(object):ENDCLASS
