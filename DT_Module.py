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

#        for istep in xrange(ctrlCI.n_timesteps):
#            for sp in pCI.names:
#                pCI.move_particles_in_uniform_fields(sp, ctrlCI)

        dt = ctrlCI.dt
        E0 = ctrlCI.E0
        for istep in xrange(ctrlCI.n_timesteps):
            # needs dt; doesn't need n_timesteps
            for sp in pCI.names:

                if self.pCI.get_species_particle_count(sp) == 0: continue # Skip if there are no particles in this species
                qmdt = self.qom[sp]*dt
                psa = pCI.pseg_arr[sp] # segmented array for this species

                psa.init_segment_loop()
                pseg = psa.get_next_segment()
                while isinstance(pseg, np_m.ndarray):
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

        for istep in xrange(ctrlCI.n_timesteps):
            # needs dt; doesn't need n_timesteps

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

#STARTCLASS
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

#STARTCLASS
class DTscratch_C(object):
    """Contains scratch-pad values
    """

    def __init__(self):
        """ Initialize variables
        """
    
        return

#class DTscratch_C(object):ENDCLASS
