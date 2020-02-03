# Particle module

"""Particle_Module treats discrete particles moving in self and external fields.
"""

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
__all__ = ['ParticleInput_C',
           'ParticleSpecies_C',
           'Particle_C',
           'ParticleMeshBoundaryConditions_C',
           'ParticleMeshSources_C',
          ]

import sys
import math
import numpy as np_m
import h5py

from Dolfin_Module import Mesh_C

#STARTCLASS
class ParticleInput_C(object):
    """Particle input class.

       Contains the variables that describe the particles. The values are
       usually set by the user in the top-level source file.

    """

    def __init__(self):

        # Usually set from ctrl.precision
        # Example: numpy.float64
        self.precision = None

        # The coordinate_system is usually the same as DTcontrol_C.coordinate_system
        self.coordinate_system = None

        # Strings that label the force components acting on the particles
        # e.g., ['x', 'y', 'z']
        self.force_components = None

        # Usually set from ctrl.precision
        # Example: numpy.float64
        self.force_precision = None

        # The flag to use either Python (default) or C++ to store and move particles
        self.use_cpp_integrators = None
        self.cpp_module = None
        
        # Values: 'loop-on-particles', 'loop-on-cells'
        self.particle_integration_loop = None

        # Determines the particle-storage dimensions
        # Example: ['x', 'y',]
        self.position_coordinates = None

# May want things like this in order to call DnT from a loop?
# or spawn off many runs?
# maybe don't need all of these:
        self.particle_species = None

        # The initial particle mesh is a copy of the field mesh
#        self.pmesh = df_m.Mesh(mesh)
#        self.pmesh_M = None

        return

#class ParticleInput_C(object):ENDCLASS

#STARTCLASS
class ParticleSpecies_C(object):
    """ParticleSpecies_C implements a particle species.

       For species-specific attributes, this can be used as the parent class.
    """

    def __init__(self, name, charge, mass, dynamics, integrator_name):
        """Initialize a ParticleSpecies_C instance.

           :cvar str name: An arbitrary name but unique for the species
           :cvar double charge: The electric charge of a single physical particle.
           :cvar double mass: The mass of a single physical particle.
           :cvar str dynamics: A string identifying how particles are pushed.  Valid
                               values are 'explicit' and 'implicit'.

           :cvar str integrator_name: A string identifying the particle-advance algorithm
                                 for this species. Valid values are
                                'advance_neutral_species'.

        """

        printWarningNonZeroCharge = True
        
        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        # Value checking
        if printWarningNonZeroCharge is True:
            if dynamics == 'neutral' and charge != 0:
                print(fncName, "\tDnT WARNING: Species %s is neutral, but has a non-zero value of electric charge: %.3g. The charge will be ignored." % (name, charge))

        self.name = name
        self.charge = charge
        self.mass = mass
        self.dynamics = dynamics
        self.integrator_name = integrator_name
        
        return
#    def __init__(self, name, charge, mass, dynamics):ENDDEF

#class ParticleSpecies_C(object):ENDCLASS

#STARTCLASS
class Particle_C(object):
    """Particle_C contains the attributes of a plasma that is represented kinetically (particle positions and velocities):
          Number of kinetic species in the plasma
          Names of the kinetic species
          Storage for the particles

       It has the attributes needed to initialize the plasma:
          List of variables associated with each particle species ('x', 'y') and their numerical types.
          Initial particle distributions specified by the user (Checks that the initial distribution functions are defined in UserParticleDistributions.py)
          Initial number of particles per cell

       :cvar position_coordinates: The labels naming the position-coordinate components.
       :vartype position_coordinates: string array
       :cvar particle_dimension: Number of spatial coordinates in a particle position.

    """

    # Static class variables

    # Particles are stored in segments of the following length:
    SEGMENT_LENGTH = 100
#    SEGMENT_LENGTH = 5
    # Bitmasks are static class variables
    DELETE_FLAG = 0b1 # the lowest bit is 1
#    TRAJECTORY_FLAG = 0b10 # the second lowest bit is 1
    TRAJECTORY_FLAG = 0b1 << 1 # the second lowest bit is 1
    # Maximum number of facet-crossing that a particle will be tracked
    # before exit() is called.
# Moved to DT_control_C    
#    MAX_FACET_CROSS_COUNT = 100
    # Initialize the unique particle ID sequence. The first particle gets ID = 0.
    # Need to parallize this variable to get unique values when using multiple
    # processors.  Maybe have a processor ID in addition?
    UNIQUE_ID_COUNTER = 0


#class Particle_C(object):
    def __init__(self, particle_input, print_flag=False):
        """Take particle data provided by the user and create the initial plasma.
        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        
        # Set local variables from passed parameters
        precision = particle_input.precision
        particleSpecies = particle_input.particle_species
        
        # Set Particle_C variables
        self.precision = precision
        self.coordinate_system = particle_input.coordinate_system
        self.use_cpp_integrators = particle_input.use_cpp_integrators
        self.force_components = particle_input.force_components
        
        # The spatial coordinates of a particle
#        self.position_coordinates = particle_input.position_coordinates
        if self.coordinate_system == 'cartesian_x' or self.coordinate_system == '1D-spherical-radius':
            self.position_coordinates = ['x',]
        elif self.coordinate_system == 'cartesian_xy':
            self.position_coordinates = ['x', 'y',]
        elif self.coordinate_system == 'cartesian_xyz':
            self.position_coordinates = ['x', 'y', 'z']
        elif self.coordinate_system == None:
            errorMsg = "Specify a coordinate system for the particles!"
            raise RuntimeError(errorMsg) # sys.exit(errorMsg)
        else:
            errorMsg = "Unknown particle coordinate system " + self.coordinate_system
            raise RuntimeError(errorMsg) # sys.exit(errorMsg)

        self.particle_dimension = len(self.position_coordinates)

        # These are the position coordinates at the start of a push
        initialPositionCoordinates = [coord+'0' for coord in self.position_coordinates]

        self.velocity_coordinates = ['u'+coord for coord in self.position_coordinates]

        phase_coordinates = self.position_coordinates + initialPositionCoordinates + self.velocity_coordinates
#        print 'particle_c... phase_coords = ', phase_coordinates

        # This is for a reference to a UserMesh_C object for particles
        self.pmesh_M = None

        # This is for a reference to a ParticleMeshBoundaryConditions_C
        # object that handles particle boundary conditions.
        self.pmesh_bcs = None

        # Count the species
        self.number_of_species = len(particleSpecies)
        if print_flag:
            print("")
            print(fncName, "\tDnT INFO: There are", self.number_of_species, " species")

        # Put the species names in a list called "species_names"
        self.species_names = [sp.name for sp in particleSpecies]

        if print_flag: print("\tDnT INFO: Species names are", self.species_names)

        # Make a lookup dictionary to get the species class from the name
        self.species_class_dict = {sp.name: sp for sp in particleSpecies}

        # Make reverse lookup dictionary species_index[] to give the species
        # index (starting from 1) given it's name.
        isp = 0
        self.species_index = {}
        for sn in self.species_names:
            isp += 1
            self.species_index[sn] = isp
            if print_flag: print(fncName, "\tDnT INFO: Species", sn, "is number", self.species_index[sn])

        # Put the user-defined plasma attributes in the following
        # dictionaries, which are indexed by the name of the species
        self.initial_distribution_type = {}
        self.initial_distribution_function = {}

        self.charge = {}
        self.mass = {}
        self.dynamics = {}
        self.integrator_names = {}
        self.qom = {}

# don't know about number_per_cell; not a fundamental number; just a particular
# initializer constant for a particular initialization

#        self.number_per_cell = {}
        self.explicit_species = []
        self.implicit_species = []
        self.neutral_species = []

        self.sap_dict = {} # This is a dictionary of SAPs. There is a key:value pair for
                           # each species.
        self.particle_count = {}

        # The names of the particle attributes and the data type of each attribute is
        # stored in the dtype dictionary "particle_dtype".
        pvars = [coord for coord in phase_coordinates]
        pvars.append('weight')

        pvartypes = [precision for var in pvars]

        # Additional particle attributes
        
        pvars.append('bitflags')
        pvartypes.append(np_m.int32)

        pvars.append('cell_index')
        pvartypes.append(np_m.int32) # The size determines how many local cells you can have.

        pvars.append('unique_ID') # ID integer starts at 0
        pvartypes.append(np_m.int32) # The size determines how many local particles
                                     # you can count.

        pvars.append('crossings') # Count of cell-crossings
        pvartypes.append(np_m.int32) # The size determines how many cell-crossings
                                     # you can count.
                                     
        particleAttributes = {'names' : pvars, 'formats': pvartypes}
        # For a particle with 3D Cartesian coordindates, the attributes stored in a
        # structure have the string names
        # ‘x’,’y’,’z’, ‘x0’,’y0’,’z0’, ‘ux’,’uy’,’uz’,
        #              ‘weight’, ‘bitflags’, ‘cell_index’, ‘unique_ID’, ‘crossings’
        self.particle_dtype = np_m.dtype(particleAttributes)

        if print_flag: print(fncName, "\tDnT INFO: Particle metadata is %s" % self.particle_dtype)

            # Make a dictionary of the dtypes
# just use ['bitflag']
#            self.bitflag_index = self.particle_dtype['names'].index('bitflags')

        for sp in particleSpecies:
            speciesName = sp.name

            # Process user input for the defining constant values of this species
            # key: 'charge'
#            if print_flag: print "(DnT INFO) sp_dict =", sp_dict
            self.charge[speciesName] = sp.charge
            if print_flag: print(fncName, "\tDnT INFO: Charge for", speciesName, "is", self.charge[speciesName])
#            if sp.charge == 0.0:
#                self.neutral_species.append(speciesName)

            # key: 'mass'
            self.mass[speciesName] = sp.mass
            self.qom[speciesName] = sp.charge/sp.mass
            if print_flag: print(fncName, "\tDnT INFO: Charge-to-mass ratio for", speciesName, "is", self.qom[speciesName])
            # key: 'number_per_cell'

# Should number_per_cell be here?  Eg if you have just some test particles
# Put it as part of the user's description of the distribution function
#            self.number_per_cell[speciesName] = sp_dict['number_per_cell']
#            if echoFlag: print 'Particle_C: number per cell for ', speciesName, ' is ', self.number_per_cell[speciesName]

            # key: 'dynamics'
            self.dynamics[speciesName] = sp.dynamics

            if sp.dynamics == 'explicit':
                self.explicit_species.append(speciesName)
            elif sp.dynamics == 'implicit':
                self.implicit_species.append(speciesName)
            elif sp.dynamics == 'neutral':
                self.neutral_species.append(speciesName)                
            else:
                errorMsg = "Unknown type of dynamics " + sp.dynamics + ' for species ' + sp_name
                raise RuntimeError(errorMsg) # sys.exit(errorMsg)

            # key: 'integrator_names'
            self.integrator_names[speciesName] = sp.integrator_name
            
            # Process user input giving the particle-variable names and types
            # for each plasma species.  Allocate initial storage
            # for particles using segmented vectors indexed by the
            # species name.
            # If using C++ SAPs:
            if self.use_cpp_integrators is True:
                # Create the SAPs in C++ (segmented_array_pair_cpp.so)
                import segmented_array_pair_solib
                import mesh_entity_arrays_solib
                # Use the C++ SegmentedArrayPair class for particle storage.  This
                # avoids having to call back to Python from C++ to manage the
                # storage.
                if self.coordinate_system == 'cartesian_x' or self.coordinate_system == '1D-spherical-radius':
                    import particle_cartesian_x_solib as cppModule
                    self.sap_dict[speciesName] = segmented_array_pair_solib.SegmentedArrayPair_cartesian_x(self.SEGMENT_LENGTH)
                elif self.coordinate_system == 'cartesian_xy':
                    import particle_cartesian_xy_solib as cppModule
                    self.sap_dict[speciesName] = segmented_array_pair_solib.SegmentedArrayPair_cartesian_xy(self.SEGMENT_LENGTH)
                elif self.coordinate_system == 'cartesian_xyz':
                    import particle_cartesian_xyz_solib as cppModule
                    self.sap_dict[speciesName] = segmented_array_pair_solib.SegmentedArrayPair_cartesian_xyz(self.SEGMENT_LENGTH)
                self.cpp_module = cppModule
            else:
            # If using Python-created SAPs:
                import SegmentedArrayPair_Module as SA_M

                # Use the Python SegmentedArrayPair class for particle storage
#                print(fncName, "\tUsing Python SAPs for particles")                
                self.sap_dict[speciesName] = SA_M.SegmentedArrayPair_C(self.SEGMENT_LENGTH, self.particle_dtype)

            # Initialize particle count for each species
            self.particle_count[speciesName] = 0

            # Initialize the unique particle ID sequence. The first particle gets ID = 0.
#            self.unique_ID_counter = 0

 #        self.user_particles_class = userParticlesClass

        ## Reference to particle initialization (add this after construction of
        ## Particle_C instance)
        self.initial_particles_dict = None

        ## Reference to particle sources (add this after construction of Particle_C
        ## instance)
        self.particle_source_dict = None

        ## Reference to particle number densities (add this after construction of
        ## Particle_C instance)
        # Set to None until the storage is allocated
#        self.dof_number_density_dict = {s: None for s in self.species_names}

        # This is for a reference to a Trajectory_C object to handle particles
        # that have the TRAJECTORY_FLAG bit turned on.  A Trajectory_C object
        # needs to be created before such particles are encountered (e.g., when
        # initial particles are created.)
        self.traj_T = None

        # This is for a reference to a ParticleHistory_C object.
        self.histories = None

        # Create a dictionary of functions associating particle-history names with
        # the functions that return the history data.
        # Could use getattr(self) here, and avoid the 'self' arg when used.
        self.history_function_dict = {'count': getattr(Particle_C,
                                                       "get_total_particle_count"),
                                      'species_count': getattr(Particle_C,
                                                               "get_species_particle_count"),
                                      'species_KE': getattr(Particle_C, "get_species_KE"),                                      
        }
        
        # An scratch ndarray for one particle is used for trajectories
        self.one_particle_arr = np_m.empty(1, dtype=self.particle_dtype)

        # A scratch array that can hold: x,y,z, (or subset)
        self.pcoord = np_m.empty(self.particle_dimension, dtype=precision)
        # A scratch array that can hold: ux,uy,uz, (or subset)
        self.pvel = np_m.empty(self.particle_dimension, dtype=precision)
        # A scratch array that can hold: x,y,z, x0,y0,z0 (or subset)
        self.pcoord2 = np_m.empty(2*self.particle_dimension, dtype=precision)
        # A scratch array for random numbers
        self.random_vals = np_m.empty(self.particle_dimension, dtype=precision)
        # A scratch array that can hold: dx,dy,dz (or subset)
        self.dx = np_m.empty(self.particle_dimension, dtype=precision)
        self.dx_in_cell = np_m.empty(self.particle_dimension, dtype=precision)
# Not scratch and not a needed class variable        self.dxFraction = 0.0 # The fraction of dx lying in the current cell
        
        # Initialize the counter for H5Part particle writes
        self.h5_step_counter = 0
        # This is used to avoid a second write on the same timestep
        self.h5_last_write_step = -1
        # A scratch buffer for H5Part. It's length may need to increase later
        self.h5_buffer_length = self.number_of_species*self.SEGMENT_LENGTH
        self.h5_buffer = np_m.empty(self.h5_buffer_length, dtype=np_m.float64)

        # Moved negE, etc., to initialize_particle_integration()
        
        # Not used yet
        self.particle_integration_loop = particle_input.particle_integration_loop

        return
#    def __init__(self, particle_input, print_flag = False):ENDDEF

#class Particle_C(object):
    def initialize_particles(self, print_flags, neg_E_field=None):
        """Generate initial particle distributions.

           This function is similar to 'add_more_particles()', which loops on
           particle source regions during the run.

        """
#    p_P = runCI.particles # abbrev for particle Class Instance

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        # Loop on initial_particles_dict

        initialParticlesDict = self.initial_particles_dict

        # Loop on the initialization methods in the dictionary

        # For 'listed' initialization, the dictionary has entries like
        # {ipName: (ipParams,)}
        # For 'function_over_region' initialization, the dictionary has entries
        # like
        # {ipName: (ipParams, ipFunc, ipRegion)}
        for ipName, ipTuple in initialParticlesDict.items():
            ipParams = ipTuple[0]
            s = ipParams['species_name']
            initialDistributionType = ipParams['initial_distribution_type']

            if print_flags[s] is True:
                print(fncName, "Initializating", ipName, "particles")
            if initialDistributionType == 'listed':
                self.create_from_list(s, print_flags[s])
            elif initialDistributionType == 'function_over_region':
                ipFunc = ipTuple[1]
                ipRegion = ipTuple[2]
                # Invoke the creation function
                time = 0.0
                ipFunc(time, ipRegion, ipParams, neg_E_field)
            else:
                errorMsg = fncName + "\tDnT ERROR: Unknown initial_distribution_type %s for species %s" % (initial_distribution_type, s)
                raise RuntimeError(errorMsg) # sys.exit(errorMsg)

        # for sp in self.species_names:
        #     initialDistributionType = self.initial_distribution_type[sp]
        #     if initialDistributionType == 'listed':
        #         self.create_from_list(sp, print_flags[sp])
        #     elif self.initial_distribution_type == 'functional':
        #         self.create_from_functions(sp, print_flags[sp])
        #     elif self.initial_distribution_type == 'particle_file':
        #         self.create_from_file(sp, print_flags[sp])
        #     else:
        #         errorMsg = "Unknown initial_distribution_type " + self.initial_distribution_type + " in Main for species " + sp
        #         sys.exit(errorMsg)

            # Compute the particle number density if storage has been allocated
#            if self.dof_number_density_dict[s] is not None:
#                self.accumulate_number_density(s)

        return
#    def initialize_particles(self, print_flags):ENDDEF

#class Particle_C(object):
    def initialize_particle_integration(self):
        """Create objects needed to advance a particle species.

           If there are charged particles, create arrays to hold the interpolated
           electric fields.

           Set the name of an integration function for each species. Either Python or
           C++ particle integration functions can be used.

        """

        className = self.__class__.__name__
        fncName = '('+__file__+') ' + className + "." + sys._getframe().f_code.co_name + '():\n'

        # If electric forces will be applied to particles, make reusable arrays for
        # computing the solved and external E at particle positions.
        # XXThese arrays can
        # be used for particle integration by both Python and C++ functionsXX.
        if self.force_components is not None:
            Ecomps = self.force_components
            force_precision = self.force_precision
            if self.use_cpp_integrators is True:
                nComps = len(Ecomps)
#                HERE
# make these have a shape, with C-ordering py::array::c_stype. pybind11 12.2.2, 12.2.5
                # "NumPy arrays are by default in row-major order".
#                self.negE = np_m.empty((self.SEGMENT_LENGTH, nComps), dtype=force_precision)
                self.negE = np_m.empty((nComps, self.SEGMENT_LENGTH), dtype=force_precision)
                # Change order of these too?
                self.Eext = np_m.empty((self.SEGMENT_LENGTH, nComps), dtype=force_precision)
                self.zeroE = np_m.empty((self.SEGMENT_LENGTH, nComps), dtype=force_precision)
                self.cpp_module.initialize_particle_integration(self.negE, self.Eext, self.zeroE)

            else:
                Etypes = [force_precision for comp in Ecomps] # variable type
                Eseg_dict = {'names': Ecomps, 'formats': Etypes}
                # reusable array for the solved field:
                self.negE = np_m.empty(self.SEGMENT_LENGTH, dtype=Eseg_dict)
                # Same thing, for an external field:
                self.Eext = np_m.empty(self.SEGMENT_LENGTH, dtype=Eseg_dict)
                # Create a reusable E array with all components zero
                self.zeroE = np_m.zeros(self.SEGMENT_LENGTH, dtype=Eseg_dict)

                # self.negE1 is a reusable array for one-particle field data for a
                # trajectory
                # Make an explicit name like 'Ex', 'Ey', 'Ez'
                E1comps = ['E'+comp for comp in Ecomps]
                E1seg_dict = {'names': E1comps, 'formats': Etypes}
                self.negE1 = np_m.empty(1, dtype=E1seg_dict)
                # Create a reusable E array with all components zero
                #            self.zeroE = np_m.zeros(len(self.negE1.dtype.fields), dtype=self.negE1.dtype[0])

                # self.Eext is a reusable array for one-particle external field data for a
                # trajectory.
                # Make an explicit name like 'Ex_ext', 'Ey_ext', 'Ez_ext'
                E1comps = ['E'+comp+'_ext' for comp in Ecomps]
                E1seg_dict = {'names': E1comps, 'formats': Etypes}
                self.Eext1 = np_m.empty(1, dtype=E1seg_dict)

        # Make a dictionary of the functions that will advance the particle species.
        self.integrators = {}
        if self.use_cpp_integrators is True:
            # Switch order?:
            for sn in self.species_names:
                # Advance the particles in this species with C++ functions.
                # Use the integrator specialized to the right number of cell-facets.
                tDim = self.pmesh_M.mesh.topology().dim()
                nFacets = tDim + 1
                integratorName = self.integrator_names[sn] + "_" + str(nFacets) + "_facets"
                integrator = getattr(self.cpp_module, integratorName)
                self.integrators[sn] = integrator
        else:
            for sn in self.species_names:
                # Advance the particles in this species with the Python functions in this module
                self.integrators[sn] = getattr(self, self.integrator_names[sn])
                
        return
#    def initialize_particle_integration(self): ENDDEF    


#class Particle_C(object):
    def create_from_list(self, species_name, print_flag=False):
        """Generates particles for a species from a list provided by the user.
        """

        printParticles = False
        printWarningNoTrajectory = False

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():'

        userParticlesClass = self.user_particles_class

        psegArrSp = self.sap_dict[species_name] # The SegmentedArrayPair_C object for this species

        # The function that has the list of particles
        p_list_method = getattr(userParticlesClass, species_name)

        # Create the particles
        number_of_macroparticles, particle_list = p_list_method(type='listed')

        # Check the length of the particle data by looking at the
        # datatype of the first segment of the segmented array.

# FIX for C++ version of psegArrSp:

##        if len(particle_list[0]) != len(psegArrSp[0].dtype):
##            errorMsg = fncName + "\tDnT ERROR: Species %s. Expecting particle data tuple %s, but data for first particle is: %s, which does not match up. Check UserParticles" %  (species_name, psegArrSp[0].dtype.names, particle_list[0])

            
#            print fncName, "Expect particle data for", psegArrSp[0].dtype.names
#            print "First particle is:", particle_list[0]
##            sys.exit(errorMsg)

        # Add the particles to storage and set up trajectories

#        print 'dtype = ', self.particle_dtype
#        print self.particle_dtype['names']
#        print 'index = ', self.particle_dtype['names'].index('bitflags')

        if printParticles or printWarningNoTrajectory:
            print(fncName)
        for i in range(number_of_macroparticles):
            if printParticles: print("\tspecies_name, particle_list[i] = ", species_name, particle_list[i])
            segIndex, fullIndex = psegArrSp.push_back(particle_list[i])
            # Check if this particle has the trajectory flag turned on
            # p = psegArrSp.get_item(fullIndex)
            # For compatibility with the C++ version of SegmentedArrays:
            (pseg, offset) = psegArrSp.get_segment_and_offset(fullIndex)
            p = pseg[offset]
            if p['bitflags'] & self.TRAJECTORY_FLAG != 0:
# or: p['bitflags'] should work?
                if self.traj_T is not None:
#                    print 'fullIndex for trajectory = ', fullIndex
                    dynamicsType = self.dynamics[species_name]
                    self.traj_T.create_trajectory(species_name, fullIndex, dynamicsType, unique_id_int=p['unique_ID'])
                else:
# Instead of printing this message, a traj_T object could be created here.
                    if printWarningNoTrajectory: print("\tDnT WARNING: A trajectory flag is on, but no trajectory object has been created yet.")

#        if (print_flag): print fncName, "weight for ", species_name, " is ", weight
#        if (print_flag): print fncName, "bitflags for ", species_name, " is ", bitflags

        return
#    def create_from_list(self, species_name, print_flag=False):ENDDEF


#class Particle_C(object):
    def create_from_functions(self, species_name, print_flag = False):
        """Generates particles for a species from a list provided by the user.
        """
#NOT FINISHED OR USED
        
        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'
        
        psegArrSp = self.sap_dict[species_name] # The SegmentedArrayPair_C object for this species

        p_list_method = self.initial_distribution_function[species_name]
        number_of_macroparticles, particle_list = p_list_method(type='listed')

        # ?check?
        number_of_real_particles = number_of_macroparticles/abs(self.charge[species_name])

        # Add the particles to storage

        for c in charge:
            new_p = p_vec[0]
            new_p['x'] = mp.x()
            new_p['y'] = mp.y()

        weight = number_of_real_particles/num_per_cell

        if (print_flag): print(fncName, "\tDnT INFO: Weight for", species_name, "is", weight)

# Lay down species in this cell, according to the specified density
# For now, put all particles at the centroid

        new_p = p_vec[0]
        new_p['x'] = mp.x()
        new_p['y'] = mp.y()

# Set the velocity
        vx = 0

# Put the values in a struct to save in the SegmentedArray:
        newpart = {'x':x, 'y':y, 'vx':vx, 'vy':vy, 'w':weight, }

# Store the particle in a SegmentedArray
        p_vec.add_item(newpart)

        return
#    def create_from_functions(self, species_name, print_flag = False):ENDDEF

#class Particle_C(object):
    def add_more_particles(self, ctrl, neg_E_field=None):
        """Add particles to the existing particle distributions.

        """

        printInfoAddMoreParticles = False
        
        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'
        
        if printInfoAddMoreParticles is True:
            print(fncName, "\tDnT INFO: Function entered at timeloop_count %d, time %.3g" % (ctrl.timeloop_count, ctrl.time))

        particleSourceDict = self.particle_source_dict

        # Loop on the source regions
        for srcName, srcTuple in particleSourceDict.items():
            # Get the physical parameters, the creation function, and the source
            # region.
            (srcParams, srcFunc, srcRegion) = srcTuple
            if ctrl.timeloop_count % srcParams['timestep_interval'] == 0:
                # Invoke the creation function
                if printInfoAddMoreParticles is True:
                    print(fncName, "\tDnT INFO: Creating new particles for source %s" % (srcName))
                srcFunc(ctrl.timeloop_count, ctrl.time, srcRegion, srcParams, neg_E_field)

        return
#    def add_more_particles(self, ctrl, print_flag=False):ENDDEF

#class Particle_C(object):
    def get_species_number(self, species_name):
        """Counts the number of physical particles for the given species.

           The number of physical particles is just the sum of the particle weights.

        """

        printSpeciesNumberFlag = False

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'
        
        # Sum up the macroparticles weights for this species
        species_number = 0.0
        sap = self.sap_dict[species_name] # segmented array for this species

        # Loop over the particle "out" segments
        (npSeg, pseg) = sap.init_out_loop()
        while isinstance(pseg, np_m.ndarray):
            species_number += np_m.sum(pseg['weight'])
            (npSeg, pseg) = sap.get_next_segment("out")

        if printSpeciesNumberFlag:
            print(fncName, "\tprintSpeciesNumberFlag=True:", species_name, "species has kinetic energy", species_number)
        
        return species_number
#    def get_species_number(self, species_name):ENDDEF


#class Particle_C(object):
    def get_species_particle_count(self, species_name, print_flag=False):
        """Counts the macroparticles for a given species.
        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        psegArrSp = self.sap_dict[species_name] # storage array for this species
        number_of_macroparticles = psegArrSp.get_number_of_items()
        if print_flag: print(fncName, "\tDnT INFO:", species_name, "has", number_of_macroparticles, "macroparticles")
# this should be done by the calling function?
#        self.particle_count[species_name] = number_of_particles

        return number_of_macroparticles
#    def get_species_particle_count(self, species_name, print_flag = False):ENDDEF

#class Particle_C(object):
    def get_total_particle_count(self, print_flag = False):
        """Counts all the macroparticles.
        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'
        psum = 0
        for sn in self.species_names:
            psegArrSp = self.sap_dict[sn] # storage array for this species
            npar = psegArrSp.get_number_of_items()
            psum += npar

        if print_flag: print(fncName, "\tDnT INFO: Total number of macroparticles in", len(self.species_names), 'species is', psum)
        return psum

#    def get_total_particle_count(self, print_flag = False):ENDDEF

#class Particle_C(object):
#    def get_species_KE(self, species_name, print_flag=False):
    def get_species_KE(self, species_name):
        """Sums the kinetic energy for a given species.
        """

        printSpeciesKEFlag = False
        
        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        # Sum up the kinetic energy of this species
        species_KE = 0.0
        sap = self.sap_dict[species_name] # segmented array for this species

        # Loop over the particle "out" segments
        (npSeg, pseg) = sap.init_out_loop()
        while isinstance(pseg, np_m.ndarray):
            # Sum of 0.5*m*u^2 for particles in this segment
            for ucomp in self.velocity_coordinates:
                species_KE += np_m.sum(pseg['weight']*pseg[ucomp]**2)
            
            (npSeg, pseg) = sap.get_next_segment("out")

        species_KE *= 0.5*self.mass[species_name]

        if printSpeciesKEFlag:
            print(fncName, "\tprintSpeciesKEFlag=True:", species_name, "species has kinetic energy", species_KE)
        
        return species_KE
#    def get_species_KE(self, species_name):ENDDEF
    

#class Particle_C(object):
    def compute_mesh_cell_indices(self):
        """Compute the cell index for each particle.

           :returns: None

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'
        
        if self.pmesh_M is not None:
            for sn in self.species_names:
# There could be a branch to a C++ particle loop here, if use_cpp_integrators is True. Then we
# wouldn't need to compute the Python version of cell_dict{} (needed by is_inside_cell()).
                sap = self.sap_dict[sn] # segmented array for this species
                (npSeg, pseg) = sap.init_out_loop()

                while isinstance(pseg, np_m.ndarray):
    #                for ip in xrange(pseg.size):
                    for ip in range(npSeg): # Could use for p in pseg[0:npSeg] instead
                        pseg[ip]['cell_index'] = self.pmesh_M.compute_cell_index(pseg[ip])
#                        print('ip, index =', ip, pseg[ip]['cell_index'])
# Check that is_inside_cell() confirms the cell index:
#                        print('pseg[ip] =', pseg[ip])
                        if not self.pmesh_M.is_inside_cell(pseg[ip], pseg[ip]['cell_index']):
#                        if not self.pmesh_M.is_inside_CPP(pseg[ip], pseg[ip]['cell_index']):
                            errorMsg = "%s\tDnT ERROR: is_inside_cell() check failed for particle %d" % (fncName, ip)
                            raise RuntimeError(errorMsg) # sys.exit(errorMsg)
#                        else:
#                            print fncName, "*** is_inside_cell check passes for particle", pseg[ip], "***"
                    (npSeg, pseg) = sap.get_next_segment("out")
        else:
            print(fncName, "\tDnT WARNING: The reference to pmesh_M is None")

        return
#    def compute_mesh_cell_indices(self):ENDDEF

#class Particle_C(object):
    def initialize_particle_mesh(self, mesh_M):
        """Compute mesh data needed to move particles.

           :param mesh_M: A Mesh_C object used by the particle movers.

           :returns: None

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'
        
        if mesh_M is not None:
            self.pmesh_M = mesh_M
            if self.use_cpp_integrators is True:
                # Create a MeshEntityArrays object, which has non-standard mesh entity
                # lists needed by particle movers.
                import mesh_entity_arrays_solib

                pmesh_df = self.pmesh_M.mesh
                # Create the name of the specialized MeshEntityArrays class with the right
                # number of cell-facets
                tDim = pmesh_df.topology().dim()
                nFacets = tDim + 1
                meaClass = "MeshEntityArrays_" + str(nFacets) + "_facets"
                meaCtor = getattr(mesh_entity_arrays_solib, meaClass)
                self.pmesh_M.mea_object = meaCtor(pmesh_df, compute_particle_mesh_maps=True)
            else:
                self.pmesh_M.compute_cell_vertices_dict() # Tested, but not used anywhere!
                self.pmesh_M.compute_cell_entity_indices_dict('facet') # Get facet indices from a cell index                
                self.pmesh_M.compute_cell_facet_normals_dict()
                self.pmesh_M.compute_cell_neighbors_dict()
                self.pmesh_M.compute_cell_volume_dict() # Compute cell volumes indexed by cell index
        else:
            errorMsg = "%s\tDnT: No particle mesh has been provided. Cannot continue." % (fncName)
            raise RuntimeError(errorMsg)
        return
#    def initialize_particle_mesh(self, mesh_M):ENDDEF

#class Particle_C(object):
    def accumulate_number_density(self, species_name, dof_number_density_F=None, cell_number_density_F=None):
        """Compute the DoF and cell number density arrays for the specified
           kinetic-particle species.

           For the DoF array, this doesn't generate an actual number-density, but rather
           the inner product vector {n*u_i*dx}, where n is the density function (a sum of
           weighted delta-function) and u_i is the shape function of the i'th DoF. The
           result needs to be divided by a volume to get a dimension of density.

           For the cell array, the particle weights are summed in each cell and divided by
           the cell volume.

           This implementation loops on particles, not on cells.

           :param species_name: The string name of a kinetic-particle species
           :param cell_number_density_F: a scalar Field_C object containing cell values.

           :returns: None

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        sap = self.sap_dict[species_name] # segmented array for this species
        (npSeg, pseg) = sap.init_out_loop()

        while isinstance(pseg, np_m.ndarray):
            for p in pseg[0:npSeg]: # Need to avoid overrunning the initialized part of the array.
                # Project the density function for this particle onto the DoF shape functions
                # to get the DoF 'density'
                if dof_number_density_F is not None:
                    dof_number_density_F.interpolate_delta_function_to_dofs(p)
                # Compute the cell number density as well, if needed.
                if cell_number_density_F is not None:
                    cell_number_density_F.add_weight_to_cell(p)
            (npSeg, pseg) = sap.get_next_segment("out")

        # Convert the cell values to a cell density    
        if cell_number_density_F is not None:
            cell_number_density_F.divide_by_cell_volumes()

        return
#    def accumulate_number_density(self, species_name):ENDDEF

#class Particle_C(object):
    def accumulate_number_density_CPP(self, species_name, dof_number_density_F=None, cell_number_density_F=None):
        """Compute the DoF and cell number-density arrays for the specified
           kinetic-particle species using C++ functions.

           For the DoF array, this doesn't generate an actual number-density, but rather
           the inner product vector {n*u_i*dx}, where n is the density function (a sum of
           weighted delta-function) and u_i is the shape function of the i'th DoF. The
           result needs to be divided by a volume to get a dimension of density.

           For the cell array, the particle weights are summed in each cell and divided by
           the cell volume.

           This implementation loops on particles, not on cells.

           :param species_name: The string name of a kinetic-particle species
           :param cell_number_density_F: a scalar Field_C object containing cell values.

           :returns: None

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        sap = self.sap_dict[species_name] # segmented array for this species
        (npSeg, pseg) = sap.init_out_loop()

        while isinstance(pseg, np_m.ndarray):
            if dof_number_density_F is not None:
            # Project the density function for the particles onto the DoF shape functions
            # to get the DoF 'density'

#  Pass pseg[0:npSeg] instead to limit the range to real particles?

                if self.particle_dimension == 3:
                    dnt_cpp.interpolate_weights_to_dofs3D(pseg, dof_number_density_F.function._cpp_object)
                elif self.particle_dimension == 2:
                    dnt_cpp.interpolate_weights_to_dofs2D(pseg, dof_number_density_F.function._cpp_object)
                elif self.particle_dimension == 1:
                    dnt_cpp.interpolate_weights_to_dofs1D(pseg, dof_number_density_F.function._cpp_object)
            # Compute the cell number density as well, if requested.
            if cell_number_density_F is not None:
                if self.particle_dimension == 3:
                    dnt_cpp.add_weights_to_cells3D(pseg, cell_number_density_F.function._cpp_object)
                elif self.particle_dimension == 2:
                    dnt_cpp.add_weights_to_cells2D(pseg, cell_number_density_F.function._cpp_object)
                elif self.particle_dimension == 1:
                    dnt_cpp.add_weights_to_cells1D(pseg, cell_number_density_F.function._cpp_object)
            (npSeg, pseg) = sap.get_next_segment("out")

        # After all the particles have been summed, convert the cell values
        # to a cell density
        if cell_number_density_F is not None:
            dnt_cpp.divide_by_cell_volumes(cell_number_density_F.function._cpp_object, self.pmesh_M.cell_volume_dict)

        return
#    def accumulate_number_density_CPP(self, species_name):ENDDEF

# Useless function:

#class Particle_C(object):
    def set_number_density(self, species_name, value):
        """Set the values in the number density array for the given species.

           :param species_name: The string name of a kinetic-particle species
           :param value: The number density is set to this value.

           :returns: None

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        sys.exit("set_number_density() should not be called")
#        self.dof_number_density_dict[species_name].set_values(value)

        return
#    def set_number_density(self, species_name, value):ENDDEF

#class Particle_C(object):
    def move_particles_in_electrostatic_field(self, ctrl, neg_E_field=None, external_E_field=None, accel_only=False):
        """Advance charged particles by one time increment in an electric field
           interpolated from a mesh.
        
           Apply interpolated electric force to particles.  Compute the
           change in velocity and position in time dt. Use an explicit
           method to integrate the orbit.

           If a particle leaves its initial cell, the cell that the particle
           moves to is calculated by finding what facet the particle crosses,
           and looking up what the neighbor cell is. This is repeated until the
           cell containing the final position is found.

           :param ctrl: A DTcontrol_C object.
           :param neg_E_field: A Field_C object containing the field -E
                               calculated by the solving the field equation.
           :param external_E_field: A Field_C object containing an external
                                    electric field.
           :param bool accel_only: If True, apply force to particles, but do not move them.

           :cvar int cFacet: The cell-level (not mesh-level) index number of the
                             facet crossed by a particle during the move.
           :cvar int pDim: Number of spatial coordinates in the particle location.

        """

        printInfoBoundaryCrossing = False

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        # Use local names for the passed parameters
        dt = ctrl.dt
        step = ctrl.timeloop_count # The number of timesteps done so far.
        time = ctrl.time # The time being integrated to.
        MAX_FACET_CROSS_COUNT = ctrl.MAX_FACET_CROSS_COUNT
        
        pmesh_M = self.pmesh_M
        pDim = self.particle_dimension

        # Scratch space
        pCoord2 = self.pcoord2 # x,y,z, x0,y0,z0 (or subset)
        dx = self.dx # dx is the distance moved in one step.

        for sn in self.explicit_species:

            # Invariant parameters
#            print 'sn = ', sn
            qmdt = self.qom[sn]*dt

#        nlocal = self.particles.get_species_particle_count(sp, echoFlag = False)
            sap = self.sap_dict[sn] # segmented array for this species

            if self.get_species_particle_count(sn) == 0: continue

            # negE is long enough to hold one particle segment worth
            # of vector field values.

#           Accelerate and move all the particles in this species
#            (npSeg, psegIn, psegOut) = sap.init_inout_loop()
#            (npSeg, psegIn) = sap.init_inout_loop()
            (npSeg, psegIn, psegOut) = sap.init_inout_loop()
            particleCount = 0
            # ipOut counts particles being written to the current
            # "out" segment.
            ipOut = 0
            indexChange = False # Flag to indicate if particle SA indices have
                                # changed.  This affects, e.g., trajectories, which
                                # use SA indices to identify particles chosen for
                                # trajectory plots.
            while isinstance(psegIn, np_m.ndarray):
#            print psegIn['z']
#                print 'particles_mod: particles module: psegIn = ', psegIn
                # Compute electric field for each particle

#                print("Before:")
#                print("position:", psegIn['x'][0:npSeg], psegIn['y'][0:npSeg]) #, psegIn['z']
#                print("velocity:", psegIn['ux'][0:npSeg], psegIn['uy'][0:npSeg]) #, psegIn['uz']

                Eseg = None # This is for the case where no fields are to be applied to the particles,
                if neg_E_field is not None:
                    if ctrl.apply_solved_electric_field is None or ctrl.apply_solved_electric_field[sn] is True:
                        # Interpolate the solved electric field to the particle positions
                        # NB: This is the negative electric field
                        neg_E_field.interpolate_field_to_points(psegIn[0:npSeg], self.negE)
                        # Truncate negE to the number of particles to get the
                        # += operations below to work: These operations need
                        # the vectors to be the same length.
                        # This is a ref, not a copy
                        Eseg = self.negE[0:npSeg]
                        # Flip the sign to get correct E. (This is a vector
                        # operation on each component of E.)
                        for n in Eseg.dtype.names:
                            Eseg[n] *= -1.0
                else:
                    self.negE[0:npSeg] = self.zeroE[0:npSeg]
                    Eseg = self.negE[0:npSeg] # A ref, not a copy.
                if external_E_field is not None:
                    if ctrl.apply_random_external_electric_field is None or ctrl.apply_random_external_electric_field[sn] is True:
                        # Interpolate the external electric field to the particle positions
                        external_E_field.interpolate_random_field_to_points(psegIn, self.Eext)
                        Eext_seg = self.Eext[0:npSeg]
                        # The following is a vector operation on each component of E. This is the best
                        # you can do for a structured array. I'm not sure that you can add arrays
                        # of E-field vectors in a single statement anyway.
                        for n in Eseg.dtype.names:
                            Eseg[n] += Eext_seg[n]

#Old way, not using in/out arrays:
                """
                # Advance the particle velocities: Apply the electric
                # forces for time dt: e.g., pseg['ux'] +=
                # qmdt*Eseg['x']


                for comp in Eseg.dtype.names:
                    ucomp = 'u'+comp
                    pseg[ucomp] += qmdt*Eseg[comp]

                # Advance the particle positions, e.g., pseg['x'] += pseg['ux']*dt
                # NON-RELATIVISTIC: v and u are the same
                for coord in self.position_coordinates:
                    coord0 = coord+'0'
                    pseg[coord0] = pseg[coord] # Save the starting positions
                    ucomp = 'u'+coord
                    pseg[coord] += pseg[ucomp]*dt

#                print "After:"
#                print "position:", pseg['x'], pseg['y'] #, pseg['z']
#                print "velocity:", pseg['ux'], pseg['uy'] #, pseg['uz']
#                print "E:", Eseg['x'], Eseg['y'] #, Eseg['z']

                """

                # Loop on the particles in this "in" segment.
                # ipIn counts through the "in" segment
                for ipIn in range(npSeg): # Use for p in pseg[0:npSeg] instead?
#                for ip in xrange(pseg.size):

                    # pseg[i] has the 'x', 'y', 'z', 'x0',..., 'ux', 'uy',... values
                    # of the ith item So pseg[i][0:3] is 'x', 'y', 'z', be we can't
                    # use slice syntax unless the array data are homogeneous.

                    # Get the next "out" segment if the current "out" segment is
                    # full. If there are no more "out" segments, allocate a new one.
                    if ipOut == self.SEGMENT_LENGTH:
                        psegOut = sap.get_next_out_segment()
                        ipOut = 0 # Reset the slot counter for the new "out" segment

                    # COPY everything from "in" to "out", to ensure weights, flags, etc. get copied.
                    psegOut[ipOut] = psegIn[ipIn] 

                    # Accelerate the particle with the field components in Eseg
                    if Eseg is not None:
                        for comp in Eseg.dtype.names:
                            ucomp = 'u'+comp
                            psegOut[ipOut][ucomp] += qmdt*Eseg[ipIn][comp] # Index OK?

# Instead of this, could do ifs on the dimension of Eseg
                    """
                    psegOut[ipOut]['ux'] = psegIn[ipIn]['ux'] + qmdt*Eseg[ipIn]['x']
                    psegOut[ipOut]['uy'] = psegIn[ipIn]['uy'] + qmdt*Eseg[ipIn]['y']
                    psegOut[ipOut]['uz'] = psegIn[ipIn]['uz'] + qmdt*Eseg[ipIn]['z']
                    """

                    # Move the particle
                    # NON-RELATIVISTIC: v and u are the same
                    if accel_only is False:
                        for coord in self.position_coordinates:
                            coord0 = coord+'0'
                            #                        psegOut[ipOut][coord0] = psegIn[ipIn][coord] # Save the starting positions
                            psegOut[ipOut][coord0] = psegOut[ipOut][coord] # Save the starting positions
                            ucomp = 'u'+coord
                            #                        psegOut[ipOut][coord] = psegIn[ipIn][coord] + psegIn[ipIn][ucomp]*dt
                            psegOut[ipOut][coord] += psegOut[ipOut][ucomp]*dt

                    """
                    psegOut[ipOut]['x0'] = psegIn[ipIn]['x']
                    psegOut[ipOut]['y0'] = psegIn[ipIn]['y']
                    psegOut[ipOut]['z0'] = psegIn[ipIn]['z']

                    psegOut[ipOut]['x'] = psegIn[ipIn]['x'] + psegOut[ipOut]['ux']*dt
                    psegOut[ipOut]['y'] = psegIn[ipIn]['y'] + psegOut[ipOut]['uy']*dt
                    psegOut[ipOut]['z'] = psegIn[ipIn]['z'] + psegOut[ipOut]['uz']*dt
                    """

                    # Check if the particle has left this cell.

#                    print 'ip, index =', ip, pseg[ip]['cell_index']
                    pCellIndex = psegOut[ipOut]['cell_index']

#                    print fncName, ": ip, pindex", ip, pCellIndex, "cell index:", pmesh_M.compute_cell_index(pseg[ip])
                    mLastFacet = Mesh_C.NO_FACET
                    facetCrossCount = 0
                    tStart = time - dt
                    dtRemaining = dt
# TODO: fix the CPP version to allow DnT_pstruct args of any dimension.
                    while not pmesh_M.is_inside_cell(psegOut[ipOut], pCellIndex):
                        # The particle has left this cell.  We
                        # need to track it across each facet in case
                        # there's a boundary-condition on that facet.
#                        print fncName, "particle has migrated"

                        facetCrossCount += 1
                        # Check for an abnormal number of facet crossings:
                        if facetCrossCount > MAX_FACET_CROSS_COUNT:
                            errorMsg = "%s\tExiting because MAX_FACET_CROSS_COUNT = %d was exceeded!" % (fncName, MAX_FACET_CROSS_COUNT)
                            raise RuntimeError(errorMsg) # sys.exit(errorMsg)

                        # Compute dx[], the move vector that starts in
                        # the current cell
                        '''
                        i=0
                        for coord in self.position_coordinates:
                            pCoord2[i] = psegOut[ipOut][coord] # Present position
                            coord0 = coord+'0'
                            pCoord2[i+pDim] = psegOut[ipOut][coord0] # Start of path in this cell
                            i+=1
                        '''
                        # Replace above by:
                        for i in range(2*pDim):
                            pCoord2[i] = psegOut[ipOut][i]

# NB: if you write "dx = pCoord2" here, then this reassigns dx, so it's not pointing
# at self.dx any more.
                        dx[:] = pCoord2[0:pDim] - pCoord2[pDim:2*pDim] # the particle move
                                                                       # vector starting in
                                                                       # the current cell.

# See save10/Particle_Module.py for a different search, where only nearby cells are looked at.

                        # Since the particle is no longer in the current cell, find which
                        # facet of the current cell it crossed.  Get (1) the cell-level index
                        # of the crossed facet, (2) the fraction of the move that's in the
                        # current cell, and (3) the vector normal to the crossed facet.
                        (cFacet, dxFraction, facetNormal) = pmesh_M.find_facet(pCoord2[pDim:2*pDim], dx, pCellIndex)
                        if cFacet != Mesh_C.NO_FACET: # If the particle crossed a facet...

                            tStart = tStart + dxFraction*dtRemaining # The starting time in the new cell
                            dtFraction = dxFraction*dtRemaining
                            dtRemaining -= dtFraction
#                            dtRemaining = (1.0 - dxFraction)*dtRemaining

#                        if found_cell != True: print "Particle is not in nearby cells; using BB search"
#                        ci = pmesh_M.compute_cell_index(pseg[ip])
#                        print "Found particle is in cell", ci
#                        pseg[ip]['cell_index'] = ci

                            # Compute the coordinates of the point where the particle
                            # crossed into the neighboring cell.  These become the
                            # new 'initial' coordinates of the particle for further
                            # tracking across facets.
#                            print fncName, "Before truncation p =", psegOut[ipOut]
                            i=0 # i indexes the move-vector coordinates (dx[i])
                            for coord in self.position_coordinates:
                                coord0 = coord+'0'
                                psegOut[ipOut][coord0] = psegOut[ipOut][coord0] + dxFraction*dx[i]
# Note: the crossing velocity could be computed here from psegIn[ipIn], and passed to record_trajectory_datum() below.                                
                                i+=1
#                            print "After truncation p =", psegOut[ipOut]

                            # Look up the mesh-level index of this facet...
                            mFacet = pmesh_M.cell_entity_indices_dict['facet'][pCellIndex][cFacet]
# mFacet should never be the same as the last facet crossed: check this
                            if mFacet == mLastFacet: # If the particle has crossed the same facet twice in succession...
                                errorMsg = "%s The mesh index of the facet crossed is %d, the same as the last facet crossed. This should not happen since the particle cannot cross the same facet twice in one move!" % (fncName, mFacet)
                                raise RuntimeError(errorMsg) # sys.exit(errorMsg)
                            else: # The particle has crossed a new facet.
                                mLastFacet = mFacet
                            # ...and get the value of the facet marker.
# mFacet is a numpy.uint32 (size_t), but the FacetFunction wants an int argument.
#                            print "type is:", type(mFacet)
                            facValue = pmesh_M.particle_boundary_marker[mFacet]
                            # Check if this facet has a non-zero marker, indicating
                            # that a callback function has to be called. E.g., the
                            # facet is a boundary.
                            if facValue != 0:
                                # Call the function associated with this value.

                                # !! A reference to dx[] exists in the BC function class.
                                if psegOut[ipOut]['bitflags'] & self.TRAJECTORY_FLAG != 0: # If this is a trajectory particle.

                                    if printInfoBoundaryCrossing is True: print("Recording boundary-crossing for particle", psegOut[ipOut]['unique_ID'])

                                    # Get the storage index that currently identifies this
                                    # particle in the trajectory list.
                                    fullIndex = sap.get_full_index(ipIn, "in")
# Why is the None here for E? Because a particle may have gone out-of-bounds?  In that case, could compute E at initial coords instead of final coords.

# Can I just send the interpolated E-field value?
                                    self.record_trajectory_datum(sn, ipOut, fullIndex, step, tStart, neg_E_field=None, external_E_field=None, facet_crossing=True)
                                self.pmesh_bcs.bc_function_dict[facValue][sn](psegOut[ipOut], sn, mFacet, dx_fraction=dxFraction, facet_normal=facetNormal)

                            # Look up the cell index of the new cell.
                            pCellIndexNew = pmesh_M.cell_neighbors_dict[pCellIndex][cFacet]
                            
                            # If the particle has left the mesh, and has been deleted, end
                            # the search.  If it hasn't been deleted (e.g., reflected),
                            # continue tracking it.
                            if pCellIndexNew == Mesh_C.NO_CELL:
                                if psegOut[ipOut]['bitflags'] & self.DELETE_FLAG == 1:
                                    psegOut[ipOut]['cell_index'] = Mesh_C.NO_CELL
                                    break # Breaks out of the facet-crossing 'while' loop.
                                # else: The boundary did not absorb the particle, so the
                                # particle is now at the cell boundary and the cell index
                                # is unchanged.
                            else:
                                psegOut[ipOut]['cell_index'] = pCellIndexNew
                                pCellIndex = pCellIndexNew # Needed for the next iteration of the while loop.
                        else: # The crossed faced is NO_FACET, which shouldn't happen.
                            errorMsg = "%s The cell index of the facet crossed is NO_FACET (%d). This should not happen since the particle has left its initial cell!" % (fncName, cFacet)
                            raise RuntimeError(errorMsg) # sys.exit(errorMsg)
#                       END:if cFacet != Mesh_C.NO_FACET:
#                   END:while not pmesh_M.is_inside_cell(psegOut[ipOut], pCellIndex)

                    # Record the number of facet-crossings
                    psegOut[ipOut]['crossings'] = facetCrossCount

                    # Don't need this since we just look up the cell
                    # index when computing negE above
                    # pseg[ip]['cell_index'] = compute_cell_index(Point(p))

                    # Check that this particle has not been deleted before
                    # incrementing the "out" particle counter.
                    if psegOut[ipOut]['bitflags'] & self.DELETE_FLAG == 0: # If this particle has not been deleted...
                        # If particle indices in the "out" SA have changed (i.e.,
                        # indexChange is True) due to deletions, then update this
                        # particle's index where needed.  E.g., if this is a
                        # trajectory particle, update its SA index in the list of
                        # trajectory particles.
                        if self.traj_T is not None:
                            if psegOut[ipOut]['bitflags'] & self.TRAJECTORY_FLAG != 0: # If this is a trajectory particle.
                                if indexChange is True:
#                                    print "mover: indexChange is true: updating the trajectory particle ID for particle", psegIn[ipIn]['unique_ID'], "ipIn =", ipIn, "ipOut =", ipOut
                                    self.update_trajectory_particleId(sn, ipIn, ipOut)

                        # Advance the "out" array counter for the next particle
                        ipOut += 1
                    else: # This particle has been deleted
                        indexChange = True # This indicates that particle SA indices
                                           # are changed past this point due to
                                           # deletions.
                        # If this was a trajectory particle, remove it's index from
                        # the trajectory-particle list.
                        if self.traj_T is not None:
                            if psegIn[ipIn]['bitflags'] & self.TRAJECTORY_FLAG != 0:
#                                print "mover: Removing particle with ID", psegIn[ipIn]['unique_ID'], "ipIn", ipIn, "from trajectories", ", ipOut is", ipOut
                                self.remove_trajectory_particleId(sn, ipIn, ipOut, step, time, dt)

                    # Check if we've reached the end of this segment.  If
                    # so, we need to start writing on a new segment.
                    if (ipOut == self.SEGMENT_LENGTH):
                        particleCount += self.SEGMENT_LENGTH

                # Done with this segment.
                # Get the next one, if it exists.
                (npSeg, psegIn) = sap.get_next_segment("in")
            # End of while isinstance(psegIn, np_m.ndarray)

            # Set values that will be used to keep track of the particle arrays
            if (ipOut != self.SEGMENT_LENGTH): # This catches the case where we exit the loop
                                               # when ipOut = SEGMENT_LENGTH and there are no
                                               # more "in" segments.  Otherwise, we would add
                                               # SEGMENT_LENGTH to particleCount twice.
                
                particleCount += ipOut
            sap.set_number_of_items("out", particleCount)
            # End of loop over segmented array

# Compute new density here?

        # End of loop over species
        return
#    def move_particles_in_electrostatic_field(self, ctrl, neg_E_field=None, external_E_field=None, accel_only=False):ENDDEF

#class Particle_C(object):
    def advance_neutral_particles(self, ctrl):
        """Advance all neutral particles by one time increment on a mesh.

           Either Python or C++ particle movers can be used.

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        for sn in self.neutral_species:
            if self.get_species_particle_count(sn) == 0: continue
            if self.use_cpp_integrators is True:
                # Advance the particles in this species with C++
                # Invoke the mover for the right number of cell-facets
                tDim = self.pmesh_M.mesh.topology().dim()
                nFacets = tDim + 1
                particleAdvancerName = "advance_neutral_species_" + str(nFacets) + "_facets"
                particleAdvancer = getattr(self.cpp_module, particleAdvancerName)
#                print(fncName, "\tCalling", particleAdvancerName)
                particleAdvancer(self, sn, ctrl)
            else:
                # Advance the particles in this species with Python
                self.advance_neutral_species(sn, ctrl)

        return
#    def advance_neutral_particles(self, ctrl):ENDDEF


#class Particle_C(object):
    def advance_neutral_species(self, species_name, ctrl, print_flag = False):
#    def integrate_neutral_species(self, species_name, ctrl, print_flag = False):
        """Advance a neutral species for one timestep using Python.

           Compute change in position in time dt. Use an explicit method to calculate the
           final position.  If a particle leaves its initial cell, the cell that the
           particle moves to is calculated by finding what facet the particle crosses, and
           looking up what the neighbor cell is. This is repeated until the cell
           containing the final position is found.

           :param ctrl: A DTcontrol_C object

           :cvar double dtRemaining: The time a particle has left to move in the new cell
                                      the particle has entered.
           :cvar int pDim: Number of spatial coordinates in the particle location.
           :cvar double tStart: The time at which a particle starts its move in the
                                 current cell.


        """

        printInfoBoundaryCrossing = False
        
        # Set local names for the passed parameters
        dt = ctrl.dt
        step = ctrl.timeloop_count
        time = ctrl.time
        MAX_FACET_CROSS_COUNT = ctrl.MAX_FACET_CROSS_COUNT

        pmesh_M = self.pmesh_M
        pDim = self.particle_dimension

        # Scratch space
        pCoord2 = self.pcoord2 # Can hold x,y,z, x0,y0,z0 (or a subset of these)
        dx = self.dx # dx is the distance moved in one step.
        
        # Invariant parameters

        qmdt = self.qom[species_name]*dt

#        nlocal = self.particles.get_species_particle_count(sp, echoFlag = False)

        # Accelerate the particles
        sap = self.sap_dict[species_name] # segmented array for this species

        # This loop moves particles, so swap the in/out particle arrays

# Also, check if the following copies, or resets the reference. Ans: it COPIES.
#                    p_out = psegOut[ipOut]
# To get a REFERENCE, need to use SLICING syntax:
#                    p_out = psegIn[ipIn][:]
# But you can't get a reference to a single element of psegIn[]. See:
# http://stackoverflow.com/questions/23654088/reference-of-a-single-numpy-array-element
# See also HPL p. 137. b=a[1,:], but b is still an array.
# What about getting a ref to one particle returned from a function? That's possible because it uses the stack?
        
        (npSeg, psegIn, psegOut) = sap.init_inout_loop() # (number of particles in
                                                         # this segment, ref to
                                                         # segment)
        particleCount = 0
        # ipOut counts particles being written to the current
        # "out" segment.
        ipOut = 0
        indexChange = False # Flag to indicate if particle SA indices have
                            # changed.  This affects, e.g., trajectories, which
                            # use SA indices to identify particles chosen for
                            # trajectory plots.
        while isinstance(psegIn, np_m.ndarray): # Keep looping until we run
                                                # out of "in" segments
            for ipIn in range(npSeg): # Loop on the particles in this "in"
                                      # segment. Could use for p in pseg[0:npSeg] instead
                # psegIn[ipIn] has the 'x', 'y', 'z', 'x0',..., 'ux', 'uy',... values of ith item
                # So psegIn[ipIn][0:3] is 'x', 'y', 'z'.
                # Can't use slice syntax here, because the array data are not of homogeneous type.

                # Skip deleted particles
                if psegIn[ipIn]['bitflags'] & self.DELETE_FLAG != 0:
                    indexChange = True # Particle SA indices are stale past
                                       # this point due to deletions.
                    continue # skip to the next particle

                # If this "out" segment is full, get the next "out" segment,
                # since there are still some "in" particles to advance.  If
                # there are no more "out" segments, allocate a new one.
                if ipOut == self.SEGMENT_LENGTH:
                    psegOut = sap.get_next_out_segment()
                    ipOut = 0 # Reset the counter for the new segment

                # The following COPIES data, rather than just setting a
                # reference to existing data.
                psegOut[ipOut] = psegIn[ipIn] # Copy this particle's data
                                              # from the input slot to the
                                              # output slot


                ## First, move the particle a full step
                # NON-RELATIVISTIC: v and u are the same
                for coord in self.position_coordinates:
                    coord0 = coord+'0'
                    psegOut[ipOut][coord0] = psegOut[ipOut][coord] # Save the starting positions
                    ucomp = 'u'+coord
                    psegOut[ipOut][coord] += psegOut[ipOut][ucomp]*dt

                #                    print "x0=", psegOut[ipOut]['x0'],"x=", psegOut[ipOut]['x']
                #                    print "psegOut", psegOut[ipOut]

                # Check if the particle has crossed into a
                # different cell. If it has, then:
                #   1. Find out which facet it crossed
                #   2. Check for a boundary-condition on that facet.
                #      a. If none (facet marker is 0), go to 3.
                #      b. If there is one, call the BC handler. If
                #      the BC allows this particle to continue,
                #      goto 3. If it doesn't, the particle is
                #      deleted.
                #
                #   3. Set the x0,y0,z0 coords to the crossing point.

                #                    print 'ip, index =', ipOut, psegOut[ipOut]['cell_index']
                pCellIndex = psegOut[ipOut]['cell_index']
                #                    print fncName, ": ip, pindex", ipOut, pCellIndex, "cell index:", pmesh_M.compute_cell_index(psegOut[ipOut])

                # Loop until the particle is in the current cell
                mLastFacet = Mesh_C.NO_FACET
                facetCrossCount = 0
                tStart = time - dt
                dtRemaining = dt

                while not pmesh_M.is_inside_cell(psegOut[ipOut], pCellIndex):
                    # The particle has left this cell.  We
                    # need to track it across each facet in case
                    # there's a boundary-condition on that facet.
                    # print (fncName, "particle has migrated")

                    facetCrossCount += 1
                    # Check for an abnormal number of facet crossings:
                    if facetCrossCount > MAX_FACET_CROSS_COUNT:
                        errorMsg = "%s !!! MAX_FACET_CROSS_COUNT exceeded!!!" % (fncName)
                        raise RuntimeError(errorMsg) # sys.exit(errorMsg)

                    # Compute dx[], the move vector that starts in
                    # the current cell
                    i=0
                    for coord in self.position_coordinates:
                        pCoord2[i] = psegOut[ipOut][coord] # Present position
                        coord0 = coord+'0'
                        pCoord2[i+pDim] = psegOut[ipOut][coord0] # Start of path in this cell
                        i+=1
                    dx = pCoord2[0:pDim] - pCoord2[pDim:2*pDim] # Move vector

                    # could return the crossing-point in pCoord2[]
                    (cFacet, dxFraction, facetNormal) = pmesh_M.find_facet(pCoord2[pDim:2*pDim], dx, pCellIndex)
                    if cFacet != Mesh_C.NO_FACET:
                    # print "facet crossed is", cFacet
                        tStart = tStart + dxFraction*dtRemaining # The starting time in the new cell
                        dtRemaining = (1.0 - dxFraction)*dtRemaining

                        # Compute the crossing point
                        # print fncName, "Before truncation p =", psegOut[ipOut]
                        i=0 # i indexes the move-vector coordinates (dx[i])
                        for coord in self.position_coordinates:
                            coord0 = coord+'0'
                            psegOut[ipOut][coord0] = psegOut[ipOut][coord0] + dxFraction*dx[i]
                            i+=1
                        # print "After truncation p =", psegOut[ipOut]

                        # Look up the mesh-level index of this facet...
                        mFacet = pmesh_M.cell_entity_indices_dict['facet'][pCellIndex][cFacet]
                        # mFacet should never be the same as the last facet crossed: check this
                        if mFacet == mLastFacet: # If the particle has crossed the same facet twice in succession...
                            errorMsg = "%s The mesh index of the facet crossed is %d, the same as the last facet crossed. This should not happen since the particle cannot cross the same facet twice in one move!" % (fncName, mFacet)
                            raise RuntimeError(errorMsg) # sys.exit(errorMsg)
                        else: # The particle has crossed a new facet.
                            mLastFacet = mFacet
                        # ...and get the value of the facet marker.
                        # mFacet is a numpy.uint32 (size_t), but the FacetFunction wants an int argument.
                        # print "type is:", type(mFacet)
                        facValue = pmesh_M.particle_boundary_marker[mFacet]
                        # Check if this facet has a non-zero marker, indicating that,
                        # e.g., the facet is a boundary.
                        if facValue != 0:
                            # Call the function associated with this value.

                            if psegOut[ipOut]['bitflags'] & self.TRAJECTORY_FLAG != 0: # If this is a trajectory particle.

                                if printInfoBoundaryCrossing is True: print("Recording boundary-crossing for particle", psegOut[ipOut]['unique_ID'])

                                # Get the storage index that currently identifies this
                                # particle in the trajectory list of particles.
# Move this call inside record_trajectory_datum:                                
                                fullIndex = sap.get_full_index(ipIn, "in")
                                self.record_trajectory_datum(species_name, ipOut, fullIndex, step, tStart, facet_crossing=True)
                            # A reference to dx[] is available in the BC function class.
                            self.pmesh_bcs.bc_function_dict[facValue][species_name](psegOut[ipOut], species_name, mFacet, dx_fraction=dxFraction, facet_normal=facetNormal)
                        # Look up the cell index of the new cell.
                        pCellIndexNew = pmesh_M.cell_neighbors_dict[pCellIndex][cFacet]

                        # If the particle has left the mesh, and has been deleted, end
                        # the search.  If it hasn't been deleted (e.g., reflected),
                        # continue tracking it.
                        if pCellIndexNew == Mesh_C.NO_CELL:
                            if psegOut[ipOut]['bitflags'] & self.DELETE_FLAG == 1:
                                psegOut[ipOut]['cell_index'] = Mesh_C.NO_CELL
                                break # Breaks out of the facet-crossing 'while' loop.
                            # else: The boundary did not absorb the particle, so the
                            # particle is now at the cell boundary and the cell index
                            # is unchanged.
                        else:
                            psegOut[ipOut]['cell_index'] = pCellIndexNew
                            pCellIndex = pCellIndexNew # Needed for the next iteration of the while loop.

                    else: # The crossed faced is NO_FACET, which shouldn't happen.
                        errorMsg = "%s The cell index of the facet crossed is %d. This should not happen since the particle has left its initial cell cell!" % (fncName, cFacet)
                        raise RuntimeError(errorMsg) # sys.exit(errorMsg)
                    # END:if cFacet != Mesh_C.NO_FACET:
                # END:while not pmesh_M.is_inside_cell(psegOut[ipOut], pCellIndex)

                # Record the number of facet-crossings
                psegOut[ipOut]['crossings'] = facetCrossCount
                # print("crossings =", facetCrossCount)

                # Check that this particle has not been deleted before
                # incrementing the "out" particle counter.
                if psegOut[ipOut]['bitflags'] & self.DELETE_FLAG == 0: # If this particle has not been deleted...
                    # If particle indices in the "out" SA have changed (i.e.,
                    # indexChange is True) due to deletions, then update this
                    # particle's index where needed.  E.g., if this is a
                    # trajectory particle, update its SA index in the list of
                    # trajectory particles.
                    if self.traj_T is not None:
                        if psegOut[ipOut]['bitflags'] & self.TRAJECTORY_FLAG != 0: # If this is a trajectory particle.
                            if indexChange is True:
                                # print "mover: indexChange is true: updating the trajectory particle ID for particle", psegIn[ipIn]['unique_ID'], "ipIn =", ipIn, "ipOut =", ipOut
                                self.update_trajectory_particleId(species_name, ipIn, ipOut)

                    # Advance the "out" array counter for the next particle
                    ipOut += 1
                else: # This particle has been deleted
                    indexChange = True # This indicates that particle SA indices
                                       # are changed past this point due to
                                       # deletions.
                    # If this was a trajectory particle, remove it's index from
                    # the trajectory-particle list.
                    if self.traj_T is not None:
                        if psegIn[ipIn]['bitflags'] & self.TRAJECTORY_FLAG != 0:
                            # print "mover: Removing particle with ID", psegIn[ipIn]['unique_ID'], "ipIn", ipIn, "from trajectories", ", ipOut is", ipOut
                            self.remove_trajectory_particleId(species_name, ipIn, ipOut, step, time, dt)

                # Check if we've reached the end of this segment.  If so, we need
                # to start writing on a new segment.  If there are no more
                # segments, allocate a new one.
                if (ipOut == self.SEGMENT_LENGTH):
                    particleCount += self.SEGMENT_LENGTH
                    # ipOut = 0 # This will cause get_next_out_segment() to be called
                               # above, if there are more "in" particles to be processed.

            # Done with this "in" segment.
            # Get the next one, if it exists.
            (npSeg, psegIn) = sap.get_next_segment("in")
        # End of while isinstance(psegIn, np_m.ndarray)                                   


        # Set values that will be used to keep track of the particle arrays
        if (ipOut != self.SEGMENT_LENGTH): # This catches the case where we exit the loop
                                           # when ipOut = SEGMENT_LENGTH and there are no
                                           # more "in" segments.  Otherwise, we would add
                                           # SEGMENT_LENGTH to particleCount twice.
            particleCount += ipOut
        sap.set_number_of_items("out", particleCount)
        # Loop over segmented array for this species ends

        return
#     def advance_neutral_species(self, species_name, ctrl, print_flag = False): ENDDEF


#class Particle_C(object):
    def move_particles_in_electrostatic_potential(self, species_names, dt, fCI, fpCI):
        """Do an energy-conserving push in the electrostatic
           potential.
        """
        for sn in species_names:

            # Invariant parameters
            qmdt = self.qom[sn]*dt

#        nlocal = self.particles.get_species_particle_count(sp, echoFlag = False)
            sap = self.sap_dict[sn] # segmented array for this species

            # Eseg is long enough to hold one particle segment worth
            # of values.  If segment length is the same for all
            # species, this could be a Field_Particle class array,
            # i.e., one that persists.

#            Eseg = np.empty(len(sap[0]), dtype=fpCI.Eseg_dict)

#           Accelerate and move all the particles in this species
            sap.init_segment_loop()
            pseg = sap.get_next_segment()
#            while pseg != None:
            while isinstance(pseg, np_m.ndarray):
#            print pseg['z']
#            print pseg
                # Compute electric field for each particle

# rewrite for potential:

#                Eseg = fpCI.compute_E_at_particles(pseg, Eseg)


                # Parabola algorithm
                # Accelerate this block of particles
#                pseg['ux'] += qmdt*Eseg['x']
#                pseg['uy'] += qmdt*Eseg['y']
#                pseg['uz'] += qmdt*Eseg['z']
                # NON-RELATIVISTIC: v and u are the same
                # Move the particles
                pseg['x'] += pseg['ux']*dt
                pseg['y'] += pseg['uy']*dt
                pseg['z'] += pseg['uz']*dt
#            print pseg['z']
                pseg = sap.get_next_segment()

        return
#    def move_particles_in_electrostatic_potential(self, species_names, dt, fCI, fpCI):ENDDEF



#class Particle_C(object):
    def move_charged_species_in_uniform_fields(self, species_name, ctrl, print_flag = False):
        """Apply the electric field ctrl.E0 to particles of the given species for a
           time interval ctrl.dt.  Compute the resulting change in particle
           velocities and positions.

           NB: There is no implementation of boundary-conditions in this function, so
           particles do not get deleted.

        """
# Use numpy array syntax to do this, instead of loops; See HPL Sec 4.2.2.

        dt = ctrl.dt
        E0 = ctrl.E0

        # Invariant parameters

        qmdt = self.qom[species_name]*dt

#        nlocal = self.particles.get_species_particle_count(sp, echoFlag = False)

        # Accelerate the particles
        sap = self.sap_dict[species_name] # segmented array for this species

        # This loop moves particles, so swap the in/out particle arrays
#        (npSeg, psegIn, psegOut) = sap.init_inout_loop()
#        (npSeg, psegIn) = sap.init_inout_loop()
        (npSeg, psegIn, psegOut) = sap.init_inout_loop()

#        print "move_charged_species_in_uniform_fields: npSeg, psegIn, psegOut:", npSeg, psegIn, psegOut
        # Get the first segment of particles to move
#        pseg = sap.get_next_segment()

        # Loop on segments until a None value is returned
        particleCount = 0
        # ipOut counts through the "out" array
        ipOut = 0
        while isinstance(psegIn, np_m.ndarray):
#            print "type pseg = ", type(pseg)
#            print pseg['z']
#            print pseg
            # Accelerate this block of particles

            # Loop on the particles in this "in" segment.
            # ipIn counts through the "in" segment
            for ipIn in range(npSeg): # Could use for p in pseg[0:npSeg] instead.

                # Get the next "out" segment if the current "out" segment is full.
                # If there are no more "out" segments, allocate a new one.
                if ipOut == self.SEGMENT_LENGTH:
                    psegOut = sap.get_next_out_segment()
                    ipOut = 0 # Reset the counter for the new segment

                    # COPY everything from "in" to "out", to ensure weights, flags, etc. get copied.
                psegOut[ipOut] = psegIn[ipIn]

                # Accelerate the particle
                psegOut[ipOut]['ux'] = psegIn[ipIn]['ux'] + qmdt*E0.x
                psegOut[ipOut]['uy'] = psegIn[ipIn]['uy'] + qmdt*E0.y
                psegOut[ipOut]['uz'] = psegIn[ipIn]['uz'] + qmdt*E0.z

                # Move the particle
                # NON-RELATIVISTIC: v and u are the same
                psegOut[ipOut]['x'] = psegIn[ipIn]['x'] + psegOut[ipOut]['ux']*dt
                # Could be written this way because of the COPY above:
                # psegOut[ipOut]['x'] += psegOut[ipOut]['ux']*dt
                psegOut[ipOut]['y'] = psegIn[ipIn]['y'] + psegOut[ipOut]['uy']*dt
                psegOut[ipOut]['z'] = psegIn[ipIn]['z'] + psegOut[ipOut]['uz']*dt

                ipOut += 1
                # Check if we've reached the end of this segment.  If
                # so, we need to start writing on a new segment.  But
                # don't increment to a new "out" segment unless there
                # are actually more "in" particles to process.
                if (ipOut == self.SEGMENT_LENGTH):
                    particleCount += self.SEGMENT_LENGTH
                    
            # Done with this segment.
            # Get the next one, if it exists.
            (npSeg, psegIn) = sap.get_next_segment("in")

        # Set values that will be used to keep track of the particle arrays
        particleCount += ipOut
        sap.set_number_of_items("out", particleCount)

        return
#    def move_charged_species_in_uniform_fields(self, species_name, ctrl, print_flag = False):ENDDEF

# Just push one segment:

    def move_particle_segment_in_uniform_fields(self, dt, qmdt, pseg, E0):
        """Push one segment of particles in fields that are uniform in space and
           time.
        """
# Use numpy array syntax to do this, instead of loops; See HPL Sec 4.2.2.

        # Accelerate this block of particles
        pseg['ux'] += qmdt*E0.x
        pseg['uy'] += qmdt*E0.y
        pseg['uz'] += qmdt*E0.z
        # Move the particles
        # NON-RELATIVISTIC: v and u are the same
        pseg['x'] += pseg['ux']*dt
        pseg['y'] += pseg['uy']*dt
        pseg['z'] += pseg['uz']*dt
#            print pseg['z']

        return
#    def move_particles_in_uniform_fields(self, species_name, ctrl, print_flag = False):ENDDEF

#class Particle_C(object):
    def update_trajectory_particleId(self, sn, ip_in, ip_out):
        """Replace the full index of a trajectory particle with the new value.

           Trajectory particles are identified by their full SA index.  When
           particles are deleted from the SA, these indices change for the
           remaining particles, and have to be updated.

           :param str sn: The name of the current species being advanced.
           :param int ip_in: The particle's index in the current "in" segment.
           :param int ip_out: The particle's index in the current "out" segment.
        """

        sap = self.sap_dict[sn] # The segmented array for this species

        # Obtain the full indices of this particle in the "in" and "out" arrays.
        (full_index_in, full_index_out) = sap.get_full_indices(ip_in, ip_out)

        # Find the position of full_index_in in the trajectory storage list.  This is
        # the full index that previously identified the particle.
        p_index = self.traj_T.particle_index_list[sn].index(full_index_in)
        # Replace this with the new full index.
        self.traj_T.particle_index_list[sn][p_index] = full_index_out

#        print "update_traj: The particle that had index", full_index_in, "now has index", full_index_out

        return full_index_in, full_index_out
#    def update_trajectory_particleId(self, sn, ip_in, ip_out):ENDDEF

#class Particle_C(object):
    def remove_trajectory_particleId(self, sn, ip_in, ip_out, step, time, dt):
        """Record the last data for a particle and remove it from the trajectory
           list.

           When a particle that is tagged as a trajectory particle is deleted,
           it has to be removed for the trajectory list.  Record the last datum
           for this particle, and then change the particle index to NO_PINDEX to
           skip it in subsequent looping over trajectory particles.

           Since this particle has been deleted, and may be out-of-bounds, the
           electric field value (if it's being recorded) is set to zero instead of
           trying to interpolate it.

           :param str sn: The name of the current species being advanced
           :param int ip_in: The particle's index in the current "in" segment
           :param int ip_out: The particle's index in the current "out" segment
           :param int step: The current simulation step
           :param float time: The current simulation time
           :param float dt: The simulation timestep

           :cvar NO_PINDEX: a flag indicating the particle no longer exists.

        """

        printWarningIndex = True
        printInfoOutOfBounds = True

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        sap = self.sap_dict[sn] # The segmented array for this species
        psegOut = sap.get_current_out_segment()
        p = psegOut[ip_out]
        
        finalStep = step
        finalTime = time

        traj_T = self.traj_T
        if finalStep == self.traj_T.last_step:
            # This particle was recorded at the beginning of this step, and has now
            # gone out of bounds after being advanced by dt, so increment the step
            # count and time for the trajectory datum.
            if printInfoOutOfBounds is True: print(fncName, "\tDnT INFO: Particle has gone out-of-bounds after the advance on timestep %d: recording its last position." % finalStep)
            finalStep += 1
            finalTime += dt
            
        p_arr = self.one_particle_arr

        # Obtain the full index of this particle. SAP keeps track of the current segment.
        full_index_in = sap.get_full_index(ip_in, "in")
#        print "remove_traj: The full index of particle", ip_in, "is", full_index_in, ". The unique_ID is", p['unique_ID']

        # Record the last position
        # Retrieve particle using its full index
#        p_arr[0] = sap.get(full_index_in)
        p_arr[0] = p

        # Copy the particle values into the trajectory

        # Find the position of full_index_in in the trajectory storage list:
        p_index = self.traj_T.particle_index_list[sn].index(full_index_in)
#tph
#        print "remove_traj: The position of particle", full_index_in, "in the traj array is", p_index
#        print "remove_traj: The current particle indices are:", self.traj_T.particle_index_list[sn]

        # Copy the particle values into the trajectory arrays for this particle
        # (The index for the new point is equal to the existing number of points
        # in the trajectory.)
        newpoint = self.traj_T.trajectory_length[sn][p_index]
        # If the trajectory array length is exceeded, return, after marking the
        # trajectory as ended.
        if newpoint > traj_T.npoints - 1:
            if printWarningIndex is True:
                print("%s\tDnT WARNING: Index %d for species %s on step %d exceeds array size %d. Point will not be recorded." % (fncName, newpoint, sn, finalStep, traj_T.npoints))
                # Mark this as a trajectory where the particle no longer exists.
            self.traj_T.particle_index_list[sn][p_index] = self.traj_T.NO_PINDEX
            return full_index_in

#        print "remove_traj: This particle is ", self.dynamics[sn]
        if self.dynamics[sn] == 'explicit':
            for comp in traj_T.explicit_dict['names']:
#                print "remove_traj: comp = ", comp
                if comp == 'step':
                    traj_T.data_list[sn][p_index][comp][newpoint] = finalStep
                elif comp == 't':
                    traj_T.data_list[sn][p_index][comp][newpoint] = finalTime
                elif comp in p_arr.dtype.names:
                    traj_T.data_list[sn][p_index][comp][newpoint] = p_arr[0][comp]
                    # The particle has been deleted and may be out-of-bounds, so set
                    # E to zero instead of interpolating.
                elif comp in self.negE1.dtype.names:
                    traj_T.data_list[sn][p_index][comp][newpoint] = 0.0
                elif comp in self.Eext1.dtype.names:
                    traj_T.data_list[sn][p_index][comp][newpoint] = 0.0
        elif self.dynamics[sn] == 'implicit':
            for comp in traj_T.implicit_dict['names']:
                print("remove_traj: comp = ", comp)
                if comp == 'step':
                    traj_T.data_list[sn][p_index][comp][newpoint] = finalStep
                elif comp == 't':
                    traj_T.data_list[sn][p_index][comp][newpoint] = finalTime
                elif comp in p_arr.dtype.names:
                    traj_T.data_list[sn][p_index][comp][newpoint] = p_arr[0][comp]
#TODO fix this for implicit species:
                # The particle has been deleted and may be out-of-bounds, so set
                # E to zero instead of interpolating.
                # elif comp in self.negE1.dtype.names:
                #     traj_T.data_list[sn][p_index][comp][newpoint] = 0.0
        elif self.dynamics[sn] == 'neutral':
            for comp in traj_T.neutral_dict['names']:
                if comp == 'step':
                    traj_T.data_list[sn][p_index][comp][newpoint] = finalStep
                elif comp == 't':
                    traj_T.data_list[sn][p_index][comp][newpoint] = finalTime
                elif comp in p_arr.dtype.names:
                    traj_T.data_list[sn][p_index][comp][newpoint] = p_arr[0][comp]
                    
        self.traj_T.trajectory_length[sn][p_index] += 1 # Increment the trajectory length

        # Mark this as a trajectory where the particle no longer exists.
        self.traj_T.particle_index_list[sn][p_index] = self.traj_T.NO_PINDEX

        return full_index_in
#    def remove_trajectory_particleId(self, sn, ip_in):ENDDEF

#class Particle_C(object):
    def record_trajectory_datum(self, species_name, ip_out, full_index, step, time, neg_E_field=None, external_E_field=None, facet_crossing=False):
        """Save a single data-record of a particle trajectory.

           NB: It is assumed that TRAJECTORY_FLAG has been checked before this
           function is called.

           :param species_name: Name of the species that the marked particle belongs
                                to.
           :param ip_out: Index of the particle in the current "out" segment.

           :param int full_index: The full storage index of the particle, which is used to identify the particle in the trajectory storage.

           :param neg_E_field: A Field_C object containing the solved field -E.

           :param external_E_field: A Field_C object containing the external
                                    electric field.

           :param bool facet_crossing: If True, this is a facet-crossing call, and the initial
                                       particle coordinates are recorded instead of the final
                                       coordinates.

        """

        printWarningIndex = True

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        # Access the particle's record in the current "out" segment. This has the
        # most up-to-date data for this particle.
        sap = self.sap_dict[species_name] # The segmented array for this species
        psegOut = sap.get_current_out_segment()
        p = psegOut[ip_out]
        # Get the storage index that currently identifies this particle in the
        # trajectory list of particles.
#        fullIndex = sap.get_full_index(ip_in, "in")
        
        traj_T = self.traj_T

        # Copy the particle into an array, since an array must be passed below.
        p_arr = self.one_particle_arr
        p_arr[0] = p

        # print "record_trajectory_datum: p_arr[0] = ", p_arr[0]

        if self.dynamics[species_name] == 'explicit':

            # Compute the SOLVED field at this single particle
            if neg_E_field is None:
                self.negE1[0] = self.zeroE[0]
            else:
# If facet_crossing is true, could interpolate to the initial coordinate, not the final.  The final may be out-of-bounds.
                
# Is there already an interpolated value available?
                neg_E_field.interpolate_field_to_points(p_arr, self.negE1)

                
            # E_arr is a ref, not a copy (HPL,p. 137)
            E_arr = self.negE1[0:p_arr.size] # Need this syntax, even though there's
                                             # just one particle, to make E_arr an array!

            # Flip the sign to get correct E. (This is a vector operation on each
            # component of E.)
            for n in E_arr.dtype.names:
                E_arr[n] *= -1.0
            # print "E_arr.dtype.names =", E_arr.dtype.names

            # Compute the EXTERNAL field at this single particle
            if external_E_field is None:
                self.Eext1[0] = self.zeroE[0]
            else:
                external_E_field.interpolate_field_to_points(p_arr, self.Eext1)

            Eext_arr = self.Eext1[0:p_arr.size] # Need this syntax, even though there's
                                                # just one particle, to make Eext_arr an array!

            # print "Eext_arr.dtype.names =", Eext_arr.dtype.names

        # Look up the trajectory index for this particle
        if traj_T is not None:
            t_idx = traj_T.particle_index_list[species_name].index(full_index)

            # Copy the particle values into the trajectory arrays for this particle
            # (The index for the new point is equal to the existing number of points
            # in the trajectory.)
            newpoint = traj_T.trajectory_length[species_name][t_idx]
            if newpoint > traj_T.npoints - 1:
                if printWarningIndex is True:
                    print("%s\tDnT WARNING: Index %d for species %s on step %d exceeds array size %d. Point will not be recorded." % (fncName, newpoint, species_name, step, traj_T.npoints))
                return

            if self.dynamics[species_name] == 'explicit':
                for comp in traj_T.explicit_dict['names']:
                    # print 'comp = ', comp
                    if comp == 'step':
                        traj_T.data_list[species_name][t_idx][comp][newpoint] = step
                    elif comp == 't':
                        traj_T.data_list[species_name][t_idx][comp][newpoint] = time
                    elif comp in p_arr.dtype.names:
                        if facet_crossing is True and comp in self.position_coordinates:
                            traj_T.data_list[species_name][t_idx][comp][newpoint] = p_arr[0][comp+'0']
                        else:
                            traj_T.data_list[species_name][t_idx][comp][newpoint] = p_arr[0][comp]
                    elif comp in E_arr.dtype.names:
                        traj_T.data_list[species_name][t_idx][comp][newpoint] = E_arr[0][comp]
                    elif comp in Eext_arr.dtype.names:
                        traj_T.data_list[species_name][t_idx][comp][newpoint] = Eext_arr[0][comp]
            elif self.dynamics[species_name] == 'neutral':
                for comp in traj_T.explicit_dict['names']:
                    # print 'comp = ', comp
                    if comp == 'step':
                        traj_T.data_list[species_name][t_idx][comp][newpoint] = step
                    elif comp == 't':
                        traj_T.data_list[species_name][t_idx][comp][newpoint] = time
                    elif comp in p_arr.dtype.names:
                        if facet_crossing is True and comp in self.position_coordinates:
                            traj_T.data_list[species_name][t_idx][comp][newpoint] = p_arr[0][comp+'0']
                        else:
                            traj_T.data_list[species_name][t_idx][comp][newpoint] = p_arr[0][comp]
            elif self.dynamics[species_name] == 'implicit':
                # TODO
                pass
                
            traj_T.trajectory_length[species_name][t_idx] += 1 # Increment the trajectory length

        return
#    def record_trajectory_datum(self, neg_E_field=None):ENDDEF

#class Particle_C(object):
    def record_trajectory_data(self, step, time, neg_E_field=None, external_E_field=None):
        """Save one trajectory data-record for all the marked particles of all
           species.

           If a particle's index in particle_index_list[] is NO_PINDEX, it means that the
           particle has been deleted from the simulation, so no further data are
           added to it's trajectory.

           :param int step: The current timestep
           :param float time: The current physical time
           :param neg_E_field: A Field_C object containing the vector field -E
                               computed from the field equation.
           :param external_E_field: A Field_C object containing the external
                                    electric field.
        """

        printInfoSecond = False
        printWarningArraySize = True
        
        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        # Avoid a second call if the timestep hasn't changed
        if step == self.traj_T.last_step:
            if printInfoSecond is True:
                print(fncName, "\tDnT INFO: Second call on timestep %d: return without storing the data." % step)
            return

        traj_T = self.traj_T

        # Scratch space
        p_arr = self.one_particle_arr
        
#        print 'p_arr.dtype.names = ', p_arr.dtype.names
#        newpoint = self.traj_T.count

#        print fncName, 'step = ', step

        # Record a trajectory data-point for particles with explicit dynamics
        for sp in traj_T.explicit_species:
            psegArrSp = self.sap_dict[sp] # The SA for this species
            for i in range(len(traj_T.particle_index_list[sp])): # i loops over the list of
                                                             # trajectory-particles for
                                                             # this species
                ip = traj_T.particle_index_list[sp][i] # Look up the full particle index
#                print "Recording data for particle", ip, "of species", sp
                
                # If a particle no longer exists, skip it.
                if ip == traj_T.NO_PINDEX: continue

                # Retrieve the particle using its full index
                # p_arr[0] = psegArrSp.get_item(ip) # pulls value from the "out" array.
                # For compatibility with the C++ version of SegmentedArrays:
                (pseg, offset) = psegArrSp.get_segment_and_offset(ip)
                p_arr[0] = pseg[offset]
                
#                print "record_trajectory_data: p_arr[0] = ", p_arr[0]

                # Compute the force on this single particle
                if neg_E_field is None:
                    self.negE1[0] = self.zeroE[0]
                else:
                    neg_E_field.interpolate_field_to_points(p_arr, self.negE1)

                E_arr = self.negE1[0:p_arr.size] # Need this syntax, even though there's
                                                 # just one particle, to make E_arr an array.

                # Flip the sign to get correct E. (This is a vector operation on each
                # component of E.)
                for n in E_arr.dtype.names:
                    E_arr[n] *= -1.0

#                print "E_arr.dtype.names =", E_arr.dtype.names

                # Compute the EXTERNAL field at this single particle
                if external_E_field is None:
                    self.Eext1[0] = self.zeroE[0]
                else:
                    external_E_field.interpolate_field_to_points(p_arr, self.Eext1)

                Eext_arr = self.Eext1[0:p_arr.size] # Need this syntax, even though there's
                                                   # just one particle, to make Eext_arr an array!
                                               
                # Copy the particle values into the trajectory
                newpoint = traj_T.trajectory_length[sp][i]
                if newpoint > traj_T.npoints - 1:
                    if printWarningArraySize is True:
                        print("%s\tDnT WARNING: Index %d for species %s on step %d exceeds array size %d. Point will not be recorded." % (fncName, newpoint, sp, step, traj_T.npoints))
                    return
                
                for comp in traj_T.explicit_dict['names']:
#                    print 'comp = ', comp
                    if comp == 'step':
                        traj_T.data_list[sp][i][comp][newpoint] = step
                    elif comp == 't':
                        traj_T.data_list[sp][i][comp][newpoint] = time
                    elif comp in p_arr.dtype.names:
                        traj_T.data_list[sp][i][comp][newpoint] = p_arr[0][comp]
                    elif comp in E_arr.dtype.names:
                        traj_T.data_list[sp][i][comp][newpoint] = E_arr[0][comp]
                    elif comp in Eext_arr.dtype.names:
                        traj_T.data_list[sp][i][comp][newpoint] = Eext_arr[0][comp]
                        
                traj_T.trajectory_length[sp][i] += 1 # Increment the trajectory length

        # Record a trajectory data-point for particles with implicit dynamics
        for sp in traj_T.implicit_species:
            psegArrSp = self.sap_dict[sp] # The SA for this species
            for i in range(len(traj_T.particle_index_list[sp])): # i loops on
                                                             # particle-index array
                ip = traj_T.particle_index_list[sp][i] # Look up the full particle index
#                print "Recording data for particle", ip, "of species", sp
                
                # If a particle no longer exists, skip
                if ip == traj_T.NO_PINDEX: continue

                # Retrieve particle using its full index
                # p_arr[0] = psegArrSp.get_item(ip) # pulls value from the "out" array.
                # For compatibility with the C++ version of SegmentedArrays:
                (pseg, offset) = psegArrSp.get_segment_and_offset(ip)
                p_arr[0] = pseg[offset]
#                print "record_trajectory_data: p_arr[0] = ", p_arr[0]

                # Compute the force on this single particle
                
# TODO: Fix for implicit species                
#                 if neg_E_field is None:
#                     self.negE1[0] = self.zeroE[0]
#                 else:
#                     neg_E_field.interpolate_field_to_points(p_arr, self.negE1)

#                 E_arr = self.negE1[0:p_arr.size] # Need this syntax, even though there's
#                                                  # just one particle, to make E_arr an array.

#                 # Flip the sign to get correct E. (This is a vector operation on each
#                 # component of E.)
#                 for n in E_arr.dtype.names:
#                     E_arr[n] *= -1.0

# #                print "E_arr.dtype.names =", E_arr.dtype.names

                # Copy the particle values into the trajectory
                newpoint = traj_T.trajectory_length[sp][i]
                for comp in traj_T.implicit_dict['names']:
#                    print 'comp = ', comp
                    if comp == 'step':
                        traj_T.data_list[sp][i][comp][newpoint] = step
                    elif comp == 't':
                        traj_T.data_list[sp][i][comp][newpoint] = time
                    elif comp in p_arr.dtype.names:
                        traj_T.data_list[sp][i][comp][newpoint] = p_arr[0][comp]
# TODO: Fix for implicit species                        
#                    elif comp in E_arr.dtype.names:
#                        traj_T.data_list[sp][i][comp][newpoint] = E_arr[0][comp]

                traj_T.trajectory_length[sp][i] += 1 # Increment the trajectory length

        # Record a trajectory data-point for particles with neutral dynamics
        for sp in traj_T.neutral_species:
            psegArrSp = self.sap_dict[sp] # The SA for this species
            for i in range(len(traj_T.particle_index_list[sp])): # i loops on
                                                             # particle-index array
                ip = traj_T.particle_index_list[sp][i] # Look up the full particle index
#                print "Recording data for particle", ip, "of species", sp
                
                # If a particle no longer exists, skip
                if ip == traj_T.NO_PINDEX: continue

                # Retrieve particle using its full index
                # For compatibility with the C++ version of SegmentedArrays:
                (pseg, offset) = psegArrSp.get_segment_and_offset(ip)
                p_arr[0] = pseg[offset]

#                print "record_trajectory_data: p_arr[0] = ", p_arr[0]

                # Copy the particle values into the trajectory
                newpoint = traj_T.trajectory_length[sp][i]
                for comp in traj_T.neutral_dict['names']:
                    if comp == 'step':
                        traj_T.data_list[sp][i][comp][newpoint] = step
                    elif comp == 't':
                        traj_T.data_list[sp][i][comp][newpoint] = time
                    elif comp in p_arr.dtype.names:
                        traj_T.data_list[sp][i][comp][newpoint] = p_arr[0][comp]
                        
                traj_T.trajectory_length[sp][i] += 1 # Increment the trajectory length
                        

        # Update "last_step" before returning
        self.traj_T.last_step = step

        return
#    def record_trajectory_data(self, neg_E_field=None):ENDDEF

#class Particle_C(object):
    def record_history_data(self, step, time):
        """Save requested particle history data for this timestep.

           The sum-over-all-species for each requested datum is calculated and stored.
           Values for each species are also calculated and stored.  
           Per-physical-particle ("ppp") values are also calculated and stored.

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

#        self.histories.data_array.dtype.names:

        counter = self.histories.counter
        
        self.histories.data_array['step'][counter] = step
        self.histories.data_array['t'][counter] = time

        # Get the values for each species separately and then sum over all the
        # species, and also compute the per-particle value.
        for hist in self.histories.scalar_histories:
            historyForAllSpecies = 0.0
            for sp in self.species_names:
                sphist = sp+'_' + hist
                # Call the function that returns the summed value.
                historyForSpecies = self.history_function_dict["species_"+hist](self, sp)
                self.histories.data_array[sphist][counter] = historyForSpecies
                # Get the number of physical particles in this species
                speciesNumber = self.get_species_number(sp)
                sp_ppp_hist = sp+'_ppp_' + hist
                # Compute and store the per-particle value
                self.histories.data_array[sp_ppp_hist][counter] = historyForSpecies/max(speciesNumber,1)
                historyForAllSpecies += historyForSpecies
            # Store to value summed over all species
            self.histories.data_array[hist][counter] = historyForAllSpecies
                
#        counter += 1 # This is not an in-place assignment.
        self.histories.counter += 1
        
        return
#    def record_history_data(self, step, time):ENDDEF

#class Particle_C(object):

# Note: This function may be superfluous. See sphere1D.py.

    def accumulate_charge_density_from_particles(self, dof_number_density_dict_F, charge_density_F):
        """Compute the charge density contributed by the kinetic particles.

           :param charge_density_F: Storage into which charge-density is accumulated.

           :returns: None

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        ## Loop over species
        if self.pmesh_M is not None:
            for s in self.species_names:
                charge = self.charge[s]
                # Loop over the number-density arrays
                charge_density_F.multiply_add(dof_number_density_dict_F[s], charge)
        else:
            errorMsg = "%s\tDnT: No particle mesh has been provided. Cannot continue." % (fncName)
            raise RuntimeError(errorMsg)

        return
#    def accumulate_charge_density_from_particles(self, dof_number_density_dict_F, charge_density_F):ENDDEF


# Standard particle distribution types

#class Particle_C(object):
    def create_maxwellian_particles(self, step, time, domain, paramDict, neg_E_field):
        """Generates particles with a Maxwellian velocty distribution and adds
           them to particle storage.

           This function loops over the cells in the given domain and generates
           particles in each one according to the parameters provided.

           Note: If a new particle is tagged as a trajectory particle, the first
           datum is recorded.

           :param time: The physical time of the simulation

           :param domain: The spatial domain in which the Maxwellian particles are added.
                          This can be the entire mesh, or a subset of it. In either case,
                          it's a list of cells.

           :param str paramDict.species_name: The name of the species getting a
                                              Maxwellian velocity distribution

           :param float paramDict.number_density: The desired number-density created in
                                                  the domain per call

           :param float paramDict.thermal_velocity: The isotropic Maxwellian thermal
                                                    velocity

           :param float[] paramDict.drift_velocity: A velocity 3-vector

           :param int paramDict.timestep_interval: Number of timesteps between calls to
                                                   this function

           :param int paramDict.number_per_cell: Number of particles per cell in the domain

        """

        printWarningNoTrajectory = True

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        cellList = domain.cell_list
        nCell = domain.ncell

        speciesName = paramDict['species_name']
        numberDensity = paramDict['number_density']
        thermalSpeed = paramDict['thermal_speed']
        driftVelocity = paramDict['drift_velocity']
        velocityCoordinateSystem = paramDict['velocity_coordinate_system']
        numberPerCell = paramDict['number_per_cell']

        # Get the storage for this species
        psegArrSp = self.sap_dict[speciesName] # The SegmentedArrayPair_C object for this species
        pDim = self.particle_dimension

        # References to scratch space
        particle = self.one_particle_arr[0] # Checked that this really is a ref.
        pCoord = self.pcoord # x, y, z (or subset)
        pVel = self.pvel # ux, uy, uz (or subset)
        random_vals = self.random_vals

        # Loop over the source domain creating particles in each cell

        pCoord[:] = 0.0 # Zero out this array because of the case where there are
                        # more particle spatial dimensions than mesh spatial
                        # dimensions

        # Multiplier for particle weight
        weightMult = numberDensity/numberPerCell
        if self.coordinate_system == '1D-spherical-radius':
            weightMult /= 4.0*np_m.pi # 4 \pi radians in a sphere

        crossings = 0 # The cell-crossing counter has to be initialized
        for icell in range(nCell):
            # Compute the particle weight (number of particles per macroparticle)
            weight = domain.volume[icell]*weightMult

            cellIndex = domain.cell_index[icell]

            # Set positions and velocities
            cellRad = domain.radius[icell]
            cellMid = domain.midpoint[icell]
            gDim = domain.gdim

            for ip in range(numberPerCell):
                # Put the first particle at the centroid
                # Put the rest at random positions in the cell
                bitflags = 0b00 # bit flags variable is all zeroes
                if ip == 0:
                    for i in range(gDim):
                        pCoord[i] = cellMid[i]
                        # Turn trajectory flag ON for the first particle in the cell:
#tph                        
#                        bitflags = bitflags | Particle_C.TRAJECTORY_FLAG
                else:
                    while True:
                        # Put the particle in a sphere
                        random_vals = np_m.random.uniform(-1.0, 1.0, gDim) # random numbers in the range [-1,1]
                        for i in range(gDim):
                            pCoord[i] = cellMid[i]+random_vals[i]*cellRad
                        # Check if its in the cell
                        if domain.in_cell(pCoord, icell):
                            break

                # Generate a velocity vector

# For 'r_theta_z' this needs to change?

                pVel[:] = np_m.random.normal(0.0, thermalSpeed, pDim) + driftVelocity

                # Fill a particle record: (x,y,z, x0,y0,z0, ux,uy,uz, weight,
                #                         bitflags, cell_index)

                # Set x[], x0[], ux[]
                for i in range(pDim):
                    particle[i] = pCoord[i]
                    particle[i+pDim] = pCoord[i]
                    particle[i+2*pDim] = pVel[i]

                # Fill in the rest of the data for this particle
                particle[3*pDim] = weight
                particle[3*pDim+1] = bitflags
                particle[3*pDim+2] = cellIndex
                particle[3*pDim+3] = Particle_C.UNIQUE_ID_COUNTER; Particle_C.UNIQUE_ID_COUNTER += 1
#tph
#                if particle[3*pDim+3] == 2448:
#                    print "create_max: adding particle with ID", particle[3*pDim+3]
#                    particle[3*pDim+1] = bitflags | Particle_C.TRAJECTORY_FLAG

                particle[3*pDim+4] = crossings
                    
                # Store the particle in the "out" segment.
                segIndex, fullIndex = psegArrSp.push_back(self.one_particle_arr)

                # If the particle is tagged as a trajectory particle, initialize its
                # trajectory information. Note: if this is being called by
                # initialize_particles(), the E-field has not been computed yet.

                if bitflags & self.TRAJECTORY_FLAG != 0:
                    if self.traj_T is not None:
                        print("A trajectory will be recorded for particle", fullIndex, "of species:", speciesName)
                        dynamicsType = self.dynamics[speciesName]
                        self.traj_T.create_trajectory(speciesName, fullIndex, dynamicsType)
                        # Record the initial datum for this new trajectory

# does this p have an ipIn, ipOut?
                        self.record_trajectory_datum(speciesName, segIndex, fullIndex, step, time, neg_E_field=None, external_E_field=None)
                    else:
# Instead of printing this message, a traj_T object could be created here?
                        if printWarningNoTrajectory is True: print(fncName, "\tDnT WARNING: A trajectory flag is on, but no trajectory object has been created yet.")

#        if (print_flag): print fncName, "weight for ", speciesName, " is ", weight
#        if (print_flag): print fncName, "bitflags for ", speciesName, " is ", bitflags

        return
#    def create_maxwellian_particles(self,):ENDDEF

#class Particle_C(object):
    def check_particle_output_parameters(self, ctrl):
        """Check the values provided by the user

           :param ctrl: A DTcontrol_C object

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        # Check the particle output attributes requested against those that have been defined

        for pA in ctrl.particle_output_attributes:
            if pA != 'species_index' and pA not in self.particle_dtype.names:
                errorMsg = fncName + "Particle attribute '" + pA + "' is not available. Available attributes are: " + str(self.particle_dtype.names)
                raise RuntimeError(errorMsg) # sys.exit(errorMsg)

        return
#    def check_particle_output_parameters(self, ctrl):ENDDEF

#class Particle_C(object):
    def initialize_particle_output_file(self, ctrl):
        """Open a H5Part files and write the header.

           :param ctrl: A DTcontrol_C object

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        self.particle_output_handle = h5File = h5py.File(ctrl.particle_output_file,"w")

        # A file is also a group: attach the following attributes
        h5File.attrs["Title:"] = ctrl.title
        h5File.attrs["Author:"] = ctrl.author

        h5File.attrs["NumParticleTypes:"] = self.number_of_species
        for s in self.species_names:
            species_index_str = str(self.species_index[s])
            # integer index of this species
            key = "particle_type_" + species_index_str
            h5File.attrs[key] = s
            # mass of this species
            key = "mass_" + species_index_str
            h5File.attrs[key] = self.mass[s]
            # charge of this species
            key = "charge_" + species_index_str
            h5File.attrs[key] = self.charge[s]

        self.h5_step_counter = 0

        # Flush the buffer
        self.particle_output_handle.flush()

        return
#    def initialize_particle_output_file(self, ctrl):ENDDEF

#class Particle_C(object):
    def write_particles_to_file(self, ctrl, close_flag=False):
        """Writes out particle attributes for one timestep.

           A new group is created when a snaphot of the particles is written.  Each
           attribute of the particles has to be in a contiguous block, so you have to
           loop through all the particles once for each attribute.

           :cvar h5Buf: Local name of the buffer used to hold the particle attributes
                        before they're written to a file

        """

        printInfoSecondCall = False
        
        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        # Avoid a second call if the timestep hasn't changed since the last call
        # (e.g., on exiting the time loop)
        step = ctrl.timeloop_count
        if step == self.h5_last_write_step:
            if printInfoSecondCall is True: print(fncName, "\tDnT INFO: Second call on timestep %d: return without writing the data." % step)
            return

        # Create a new group for this timestep
        groupName = "Step#" + str(self.h5_step_counter)
        group = self.particle_output_handle.create_group(groupName)

        # Write the simulation time
        group.attrs["TimeValue"] = ctrl.time

        ### Count the total number of particles and see if h5_buffer_length needs to
        ### be increased
        totalParticleCount = self.get_total_particle_count(print_flag = False)
#        print fncName, "totalParticleCount =", totalParticleCount, "h5_buffer_length =", self.h5_buffer_length
        if totalParticleCount > self.h5_buffer_length:
            newNumSegs = int(1 + totalParticleCount/self.SEGMENT_LENGTH)
#            print "newNumSegs =", newNumSegs
            # Allocate a new buffer
            self.h5_buffer_length = newNumSegs*self.SEGMENT_LENGTH
            self.h5_buffer = np_m.empty(self.h5_buffer_length, dtype=np_m.float64)
        h5Buf = self.h5_buffer

        ### Loop on the particle attributes to be written out
        for pA in ctrl.particle_output_attributes:

            # Make the dtype of h5Buf match the dtype of the attribute
            # H5Part seems to be all 64-bit
            if pA == 'species_index':
                h5Buf.dtype = np_m.int64
            else:
# Changed because of FutureWarning:                
#                if np_m.issubdtype(self.particle_dtype[pA], np_m.float):
                if np_m.issubdtype(self.particle_dtype[pA], np_m.float64):
                    h5Buf.dtype = np_m.float64
                elif np_m.issubdtype(self.particle_dtype[pA], np_m.integer):
                    h5Buf.dtype = np_m.int64
                else:
                    errorMsg = fncName + "The type of particle attribute " + pA + " is not float or integer. It is " + str(self.particle_dtype[pA])
                    raise RuntimeError(errorMsg) # sys.exit(errorMsg)                    

            aOff = 0 # offset into h5Buf
            ## Loop on the species
            for s in self.species_names:
                # Skip to next species if there are no particles of this species
                if self.get_species_particle_count(s) == 0: continue

                if pA == 'species_index':
                    # The 'species_index' is the same for every particle in this
                    # species, so load the value into the buffer in one go.
                    speciesIndex = self.species_index[s]
                    npSpecies = self.get_species_particle_count(s, print_flag = False)
                    h5Buf[aOff:aOff+npSpecies] = speciesIndex
#                    print "Range of h5Buf is:", aOff, aOff+npSpecies, "shape =", h5Buf.shape
#                    print "h5Buf is:", h5Buf[aOff:aOff+npSpecies]
                    aOff += npSpecies
                else:
                    ## Loop over the segments of the "out" array to get the attribute pA
                    sap = self.sap_dict[s] # segmented array for this species
                    (npSeg, pseg) = sap.init_out_loop()
#                    particleCount = 0
                    while isinstance(pseg, np_m.ndarray):
#                        print "Array for", pA, "is: ", pseg[pA], "shape = ", pseg.shape
#                        print "Range of h5Buf is:", aOff, aOff+npSeg, "shape =", h5Buf.shape
                        h5Buf[aOff:aOff+npSeg] = pseg[pA] # Copy all the values of attribute pA to
                                                          # contiguous locations in the buffer.
                        aOff += npSeg # Advance the offset by the number of values just copied.
                        (npSeg, pseg) = sap.get_next_segment("out")

#            print "h5Buf is:", h5Buf[0:aOff]
            dset = group.create_dataset(pA, data=h5Buf[0:aOff])
#                    group.create_dataset(pA, data=h5Buf[:])
#                    d2 = group.create_dataset("y",data=y[:], dtype='f')

        self.h5_step_counter += 1

        # Update "last_write_step" before returning
        self.h5_last_write_step = step

        # Flush the buffer
        self.particle_output_handle.flush()

        if close_flag is True:
            self.particle_output_handle.close()

        return
#    def write_particles_to_file(self):ENDDEF

#class Particle_C(object):ENDCLASS

#STARTCLASS
class ParticleMeshBoundaryConditions_C(object):
    """ParticleMeshBoundaryConditions_C sets up a 2D dictionary
       (bc_function_dict, indexed by boundary and species) to treat
       kinetic particles incident on a mesh boundary.

       The functions themselves are provided by the user in a
       UserParticleBoundaryFunctions_C object.
    """

    # Static class variables

    # Particle boundary-conditions are labeled by non-zero bits:
    ABSORB  = 0b1
    REFLECT = 0b1 << 1
    NUMBER_OF_STANDARD_BCs = 2

#    ISEE    = 0b1 << 2        # Ion-stimulated electron emission
#    SEE     = 0b1 << 3        # Secondary-electron emission


#class ParticleMeshBoundaryConditions_C(object):
# Look for specific boundary conditions
    def __init__(self, speciesNames, pmesh_M, userParticlesMeshFunctionsClass, print_flag = False):
        """Initialize particle callback functions (boundary conditions).

           The following function naming-scheme is used:

           default_bc(): The global default function.
           default_bc_at_name(): The default called for all species at
                                 the boundary 'name'.
           bc_at_name_for_species: The function called for 'species'
                                   crossing 'name'.

           The most specific function found for a given boundary and
           species is used.

        """
        # Set local names from passed parameters
        pBDict = pmesh_M.particle_boundary_dict

        # We want to associate the names of the user-supplied boundary functions
        # with the int values of the facet tags.  First, we swap the BC
        # dictionary keys and values: e.g., 'xmin': 1 to 1: 'xmin'. This lets us
        # use the value of the facet function, a uint (size_t), as an index to
        # get the string name of the boundary.

        pBDictInv = {v: k for k, v in pBDict.items()}

        particleBoundaryTags = list(pBDictInv.keys())

        # Initialize a new dictionary that will contain the name of a boundary
        # function for each [boundary tag][species] pair (i.e., the key is a 2D object).
        self.bc_function_dict = dict((intTag, dict((sp, None) for sp in speciesNames)) for intTag in particleBoundaryTags)

        # Find the global default boundary function (if there is one),
        # in the supplied UserParticleBoundaryFunctions object.
        bcFunctionName = 'default_bc'
        if hasattr(userParticlesMeshFunctionsClass, bcFunctionName):
            bcGlobalDefaultFunction = getattr(userParticlesMeshFunctionsClass, bcFunctionName)
        else:
            bcGlobalDefaultFunction = None

        # Loop on particle boundary tags and on particle species
        # names to find the most specific BC.

        # Overwrite the global default function with a function specific to
        # each boundary, if there is one.
        for intTag in particleBoundaryTags:
            bcFunctionName = 'default_bc_at_' + pBDictInv[intTag]
            if hasattr(userParticlesMeshFunctionsClass, bcFunctionName):
                bcBoundaryDefaultFunction = getattr(userParticlesMeshFunctionsClass, bcFunctionName)
            else:
                bcBoundaryDefaultFunction = bcGlobalDefaultFunction
            for sp in speciesNames:
                # Overwrite the default for this boundary with a
                # function specific to this boundary and species, if
                # there is one.
                bcFunctionName = 'bc_at_' + pBDictInv[intTag] + '_for_' + sp
                if hasattr(userParticlesMeshFunctionsClass, bcFunctionName):
                    bcFunction = getattr(userParticlesMeshFunctionsClass, bcFunctionName)
                else:
                    bcFunction = bcBoundaryDefaultFunction
                if bcFunction is None:
                    print("ParticleMeshBoundaryConditions_C: No boundary condition specified for", pBDictInv[intTag], "/", sp)
                elif print_flag:
                    print("ParticleMeshBoundaryConditions_C: Boundary condition for", pBDictInv[intTag], "/", sp, "is", bcFunction)
                self.bc_function_dict[intTag][sp] = bcFunction
        return
#    def __init__(self, particle_input, particle_P, print_flag = False):ENDDEF

    def absorb(self):
        return

#    def absorb(self):ENDDEF

    def reflect(self):
        return

#    def reflect(self):ENDDEF

#class ParticleMeshBoundaryConditions_C(object):ENDCLASS

#STARTCLASS
class ParticleMeshSources_C(object):
    """ParticleMeshSources_C sets up a 2D dictionary
       (srcFunctionDict, indexed by source name and species) to create
       kinetic particles on a subdomain of the particle mesh.

       The functions themselves are provided by the user in a
       UserParticleMeshSourceFunctions_C object.
    """

    # Static class variables

#class ParticleMeshSources_C(object):
    def __init__(self, speciesNames, pmesh_M, userParticlesMeshSourceFunctionsClass, print_flag = False):
        """Initialize particle source functions.

        """
        # Set local names from passed parameters
        pSrcDict = pmesh_M.particle_source_dict

        # We want to associate the names of the user-supplied cell functions with the int
        # values of the cell tags.  First, we swap the source dictionary keys and values:
        # e.g., 'xmin': 1 to 1: 'xmin'. This lets us use the cell-tag as an index to get
        # the name of the source.

        pSrcDictInv = {v: k for k, v in pSrcDict.items()}

        particleSourceTags = list(pSrcDictInv.keys())

        # Initialize a new dictionary that will contain the name of a
        # boundary function for each (boundary, species) pair.
        # The first index is the integer tag that
        # marks the boundary facets.
        self.srcFunctionDict = dict((intTag, dict((sp, None) for sp in speciesNames)) for intTag in particleSourceTags)

        # Find the global default boundary function, if there is one,
        # in the supplied UserParticleBoundaryFunctions object.
        srcFunctionName = 'default_bc'
        if hasattr(userParticlesMeshSourceFunctionsClass, srcFunctionName):
            bcGlobalDefaultFunction = getattr(userParticlesMeshSourceFunctionsClass, srcFunctionName)
        else:
            bcGlobalDefaultFunction = None

        # Loop on particle boundary tags and on particle species
        # names to find the most specific BC.

        # Overwrite the global default function with a function specific to
        # each boundary, if there is one.
        for intTag in particleSourceTags:
            srcFunctionName = 'default_bc_at_' + pSrcDictInv[intTag]
            if hasattr(userParticlesMeshSourceFunctionsClass, srcFunctionName):
                bcBoundaryDefaultFunction = getattr(userParticlesMeshSourceFunctionsClass, srcFunctionName)
            else:
                bcBoundaryDefaultFunction = bcGlobalDefaultFunction
            for sp in speciesNames:
                # Overwrite the default for this boundary with a
                # function specific to this boundary and species, if
                # there is one.
                srcFunctionName = 'bc_at_' + pSrcDictInv[intTag] + '_for_' + sp
                if hasattr(userParticlesMeshSourceFunctionsClass, srcFunctionName):
                    srcFunction = getattr(userParticlesMeshSourceFunctionsClass, srcFunctionName)
                else:
                    srcFunction = bcBoundaryDefaultFunction
                if srcFunction is None:
                    print("ParticleMeshSources_C: No boundary condition specified for", pSrcDictInv[intTag], "/", sp)
                elif printTag:
                    print("ParticleMeshSources_C: Boundary condition for", pSrcDictInv[intTag], "/", sp, "is", srcFunction)
                self.srcFunctionDict[intTag][sp] = srcFunction

        return
#    def __init__(self, particle_input, particle_P, print_flag = False):ENDDEF

#class ParticleMeshSources_C(object):ENDCLASS
