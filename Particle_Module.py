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
#        self.pmesh = df_M.Mesh(mesh)
#        self.pmesh_M = None

        # Module containing user-supplied particle distributions and
        # boundary-conditions.
        self.user_particles_module = None
        # The class containing distribution functions
        self.user_particles_class = None
        # The class containing particle boundary conditions
        self.user_particle_bcs_class = None

        return

#class ParticleInput_C(object):ENDCLASS

class ParticleSpecies_C(object):
    """ParticleSpecies_C implements a particle species.

       For species-specific attributes, this can be used as the parent class.
    """

    def __init__(self, name, charge, mass, dynamics):
        """Initialize a ParticleSpecies_C instance.

           :cvar str name: An arbitrary name but unique for the species
           :cvar double charge: The electric charge of a single physical particle.
           :cvar double mass: The mass of a single physical particle.
           :cvar str dynamics: A string identifying how particles are pushed.  Valid
                               values are 'explicit' and 'implicit'.

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'
        self.name = name
        self.charge = charge
        self.mass = mass
        self.dynamics = dynamics

        return

#class ParticleSpecies_C(object):ENDCLASS

import SegmentedArrayPair_Module as SA_M

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
    # Maximum number of facet-crossing that a particle will be tracke
    # before exit() is called.
    MAX_FACET_CROSS_COUNT = 100


#class Particle_C(object):
    def __init__(self, particleInputCI, print_flag = False):
        """Take a list of kinetic species provided by the user and create the initial plasma
        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        # Set local variables from passed parameters
        precision = particleInputCI.precision
        self.precision = precision
        self.coordinate_system = particleInputCI.coordinate_system
        particleSpecies = particleInputCI.particle_species

        # The spatial coordinates of a particle
        self.position_coordinates = particleInputCI.position_coordinates
        self.particle_dimension = len(self.position_coordinates)

        # These are the position coordinates at the start of a push
        initialPositionCoordinates = [coord+'0' for coord in self.position_coordinates]

        velocityCoordinates = ['u'+coord for coord in self.position_coordinates]

        phase_coordinates = self.position_coordinates + initialPositionCoordinates + velocityCoordinates
#        print 'particle_c... phase_coords = ', phase_coordinates

        # This is for a reference to a UserMesh_C object for particles
        self.pmesh_M = None

        # This is for a reference to a ParticleMeshBoundaryConditions_C
        # object that handles particle boundary conditions.
        self.pmesh_bcCI = None

        # Count the species
        self.number_of_species = len(particleSpecies)
        if print_flag:
            print ""
            print fncName, "(DnT INFO) There are:", self.number_of_species, " species"

        # Put the species names in a list called "species_names"
        self.species_names = [sp.name for sp in particleSpecies]

        if print_flag: print "(DnT INFO) Species names:", self.species_names

        # Make a lookup dictionary to get the species class from the name
        self.species_class_dict = {sp.name: sp for sp in particleSpecies}

        # Make reverse lookup dictionary species_index[] to give the species
        # index (starting from 1) given it's name.
        isp = 0
        self.species_index = {}
        for sn in self.species_names:
            isp += 1
            self.species_index[sn] = isp
            if print_flag: print "(DnT INFO) Species", sn, "is number", self.species_index[sn]

        # Put the user-defined plasma attributes in the following
        # dictionaries, which are indexed by the name of the species
        self.initial_distribution_type = {}
        self.initial_distribution_function = {}

        self.charge = {}
        self.mass = {}
        self.dynamics = {}
        self.qom = {}

# don't know about number_per_cell; not a fundamental number; just a particular
# initializer constant for a particular initialization

#        self.number_per_cell = {}
        self.explicit_species = []
        self.implicit_species = []
        self.neutral_species = []

        self.pseg_arr = {}
        self.particle_count = {}

        # The names of the particle attributes and the data type of each attribute is
        # stored in the dtype dictionary "particle_dtype".
        pvars = [coord for coord in phase_coordinates]
        pvars.append('weight')

        pvartypes = [precision for var in pvars]

        pvars.append('bitflags')
        pvartypes.append(np_m.int32)

        pvars.append('cell_index')
        pvartypes.append(np_m.int32) # The size determines how many local cells you can have.

#        self.particle_dtype = {'names' : pvars, 'formats': pvartypes}
        particleAttributes = {'names' : pvars, 'formats': pvartypes}
        self.particle_dtype = np_m.dtype(particleAttributes)

        if print_flag: print "(DnT INFO) Particle metadata = %s" % self.particle_dtype

            # Make a dictionary of the dtypes
# just use ['bitflag']
#            self.bitflag_index = self.particle_dtype['names'].index('bitflags')

        for sp in particleSpecies:
            speciesName = sp.name

            # Process user input for the defining constant values of this species
            # key: 'charge'
#            if print_flag: print "(DnT INFO) sp_dict =", sp_dict
            self.charge[speciesName] = sp.charge
            if print_flag: print "(DnT INFO) Charge for", speciesName, "is", self.charge[speciesName]
            if sp.charge == 0.0:
                self.neutral_species.append(speciesName)

            # key: 'mass'
            self.mass[speciesName] = sp.mass
            self.qom[speciesName] = sp.charge/sp.mass
            if print_flag: print "(DnT INFO) Charge-to-mass ratio for", speciesName, "is", self.qom[speciesName]
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
            else:
                errorMsg = "Unknown type of dynamics " + sp.dynamics + ' for species ' + sp_name

            # Process user input giving the particle-variable names and types
            # for each plasma species.  Allocate initial storage
            # for particles using segmented vectors indexed by the
            # species name.

#            self.pseg_arr[speciesName] = SA_M.SegmentedArray_C(segment_length, metadata)
            self.pseg_arr[speciesName] = SA_M.SegmentedArray_C(self.SEGMENT_LENGTH, self.particle_dtype)

            # Initialize particle count for each species
            self.particle_count[speciesName] = 0

#        self.user_particle_class = userParticleClass

        ## Reference to particle initialization (add this after construction of
        ## Particle_C instance)
        self.initial_particles_dict = None

        ## Reference to particle sources (add this after construction of Particle_C
        ## instance)
        self.particle_source_dict = None

        ## Reference to particle number densities (add this after construction of
        ## Particle_C instance)
        # Set to None until the storage is allocated
        self.number_density_dict = {s: None for s in self.species_names}

        # This is for a reference to a Trajectory_C object to handle particles
        # that have the TRAJECTORY_FLAG bit turned on.  A Trajectory_C object
        # needs to be created before such particles are encountered (e.g., when
        # initial particles are created.)
        self.traj_T = None

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

        # Initialize the counter for H5Part particle writes
        self.h5_step_counter = 0
        # A scratch buffer for H5Part. It's length may need to increase later
        self.h5_buffer_length = self.number_of_species*self.SEGMENT_LENGTH
        self.h5_buffer = np_m.empty(self.h5_buffer_length, dtype=np_m.float64)

        # Make a reusable array "self.negE" for computing -E at particle positions
        if particleInputCI.force_components is not None:
            Ecomps = particleInputCI.force_components
            force_precision = particleInputCI.force_precision
            Etypes = [force_precision for comp in Ecomps]
            Eseg_dict = {'names': Ecomps, 'formats': Etypes}
            self.negE = np_m.empty(self.SEGMENT_LENGTH, dtype=Eseg_dict)
            # self.negE1 is used for one-particle field data for a trajectory
            # Make an explicit name like 'Ex', 'Ey', 'Ez'
            E1comps = ['E'+comp for comp in Ecomps]
            E1seg_dict = {'names': E1comps, 'formats': Etypes}
            self.negE1 = np_m.empty(1, dtype=E1seg_dict)
            # Create an E with all components zero
            self.zeroE = np_m.zeros(len(self.negE1.dtype.fields), dtype=self.negE1.dtype[0])

        # Not used yet
        self.particle_integration_loop = particleInputCI.particle_integration_loop

        return

#class Particle_C(object):
    def initialize_particles(self, print_flags, neg_E_field=None):
        """Generate initial particle distributions.

           This function is similar to 'add_more_particles()', which loops on
           particle source regions during the run.

        """
#    pCI = runCI.particles # abbrev for particle Class Instance

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        # Loop on initial_particles_dict

        initialParticlesDict = self.initial_particles_dict

        # Loop on the initialization methods in the dictionary

        # For 'listed' initialization, the dictionary has entries like
        # {ipName: (ipParams,)}
        # For 'function_over_region' initialization, the dictionary has entries
        # like
        # {ipName: (ipParams, ipFunc, ipRegion)}
        for ipName, ipTuple in initialParticlesDict.iteritems():
            ipParams = ipTuple[0]
            s = ipParams['species_name']
            initialDistributionType = ipParams['initial_distribution_type']

            if print_flags[s] is True:
                print fncName, "Initializating", ipName, "particles"
            if initialDistributionType == 'listed':
                self.create_from_list(s, print_flags[s])
            elif initialDistributionType == 'function_over_region':
                ipFunc = ipTuple[1]
                ipRegion = ipTuple[2]
                # Invoke the creation function
                time = 0.0
                ipFunc(time, ipRegion, ipParams, neg_E_field)
            else:
                errorMsg = "(DnT ERROR) %s Unknown initial_distribution_type %s for species %s" % (fncName, initial_distribution_type, s)
                sys.exit(errorMsg)

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
            if self.number_density_dict[s] is not None:
                self.accumulate_number_density(s)

        return
#    def initialize_particles(self, print_flags):ENDDEF

#class Particle_C(object):
    def create_from_list(self, species_name, print_flag=False):
        """Generates particles for a species from a list provided by the user.
        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'
        userParticleClass = self.user_particle_class

        pseg_arrCI = self.pseg_arr[species_name] # The SegmentedArray_C object for this species

        # The function that has the list of particles
        p_list_method = getattr(userParticleClass, species_name)

        # Create the particles
        number_of_macroparticles, particle_list = p_list_method(type='listed')

        # Check the length of the particle data by looking at the
        # datatype of the first segment of the segmented array.

        if len(particle_list[0]) != len(pseg_arrCI[0].dtype):
            errorMsg = "(DnT ERROR) %s Species %s. Expecting particle data tuple %s, but data for first particle is: %s, which does not match up. Check UserParticles" %  (fncName, species_name, pseg_arrCI[0].dtype.names, particle_list[0])
#            print fncName, "Expect particle data for", pseg_arrCI[0].dtype.names
#            print "First particle is:", particle_list[0]
            sys.exit(errorMsg)

        # Add the particles to storage and set up trajectories

#        print 'dtype = ', self.particle_dtype
#        print self.particle_dtype['names']
#        print 'index = ', self.particle_dtype['names'].index('bitflags')

        for i in range(number_of_macroparticles):
#            print 'species_name, particle_list[i] = ', species_name, particle_list[i]
            p, pindex = pseg_arrCI.put(particle_list[i])
            # Check if this particle has the trajectory flag turned on
            if p['bitflags'] & self.TRAJECTORY_FLAG != 0:
# or: p['bitflags'] should work?
                if self.traj_T is not None:
#                    print 'pindex for trajectory = ', pindex
                    dynamicsType = self.dynamics[species_name]
                    self.traj_T.create_trajectory(species_name, pindex, dynamicsType)
                else:
# Instead of printing this message, a traj_T object could be created here.
                    print fncName, "(DnT WARNING) A trajectory flag is on, but no trajectory object has been created yet."

#        if (print_flag): print fncName, "weight for ", species_name, " is ", weight
#        if (print_flag): print fncName, "bitflags for ", species_name, " is ", bitflags

        return
#    def create_from_list(self, species_name, print_flag=False):ENDDEF


#class Particle_C(object):
    def create_from_functions(self, species_name, print_flag = False):
        """Generates particles for a species from a list provided by the user.
        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'
        pseg_arrCI = self.pseg_arr[species_name] # The SegmentedArray_C object for this species

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

        if (print_flag): print fncName, "(DnT INFO) Weight for", species_name, "is", weight

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
    def add_more_particles(self, ctrl, neg_E_field=None, print_flag=False):
        """Add particles to the existing particle distributions.

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'
        if print_flag is True:
            print "(DnT INFO) %s\tFunction entered at timestep_count %d, time %.3g" % (fncName, ctrl.timestep_count, ctrl.time)

        particleSourceDict = self.particle_source_dict

        # Loop on the source regions
        for srcName, srcTuple in particleSourceDict.iteritems():
            # Get the physical parameters, the creation function, and the source
            # region.
            (srcParams, srcFunc, srcRegion) = srcTuple
            if ctrl.timestep_count % srcParams['timestep_interval'] == 0:
                # Invoke the creation function
                if print_flag is True:
                    print "(DnT INFO) %s\t Creating new particles for source %s" % (fncName, srcName)
                srcFunc(ctrl.timestep_count, ctrl.time, srcRegion, srcParams, neg_E_field)

        return
#    def add_more_particles(self, ctrl, print_flag=False):ENDDEF

#class Particle_C(object):
    def get_species_particle_count(self, species_name, print_flag=False):
        """Counts the particles for a given species.
        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        pseg_arrCI = self.pseg_arr[species_name] # storage array for this species
        number_of_particles = pseg_arrCI.get_number_of_items()
        if print_flag: print fncName, "(DnT INFO)", species_name, "has", number_of_particles, "macroparticles"
# this should be done by the calling function?
#        self.particle_count[species_name] = number_of_particles

        return number_of_particles
#    def get_species_particle_count(self, species_name, print_flag = False):ENDDEF

#class Particle_C(object):
    def get_total_particle_count(self, print_flag = False):
        """Counts all the particles.
        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'
        psum = 0
        for sn in self.species_names:
            pseg_arrCI = self.pseg_arr[sn] # storage array for this species
            npar = pseg_arrCI.get_number_of_items()
            psum += npar

        if print_flag: print fncName, ", (DnT INFO) Total number of macroparticles in", len(self.species_names), 'species is', psum
        return psum

#    def get_total_particle_count(self, print_flag = False):ENDDEF

#class Particle_C(object):
    def compute_mesh_cell_indices(self):
        """Compute the cell index for each particle.

           :returns: None

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'
        if self.pmesh_M is not None:
            for sn in self.species_names:
                psaCI = self.pseg_arr[sn] # segmented array for this species
                (npSeg, pseg) = psaCI.init_out_loop()

                while isinstance(pseg, np_m.ndarray):
    #                for ip in xrange(pseg.size):
                    for ip in xrange(npSeg):
                        pseg[ip]['cell_index'] = self.pmesh_M.compute_cell_index(pseg[ip])
    #                    print 'ip, index =', ip, pseg[ip]['cell_index']
# Check that is_inside() confirms the cell index:
                        if not self.pmesh_M.is_inside(pseg[ip], pseg[ip]['cell_index']):
                            errorMsg = "%s (DnT ERROR) is_inside() check failed for particle %d" % (fncName, ip)
                            sys.exit(errorMsg)
#                        else:
#                            print fncName, "*** is_inside check passes for particle", pseg[ip], "***"
                    (npSeg, pseg) = psaCI.get_next_segment('out')
        else:
            print fncName, "(DnT WARNING) The reference to pmesh_M is None"

        return
#    def compute_mesh_cell_indices(self):ENDDEF

#class Particle_C(object):
    def accumulate_number_density(self, species_name):
        """Add contributions to the number density for the specified kinetic-particle
           species.

           This implementation doesn't generate an actual number-density, but
           rather the vector {n*u_i*dx}, where the {u_i} are the shape functions
           for the i'th DoF. The result needs to be divided by a volume to get an
           actual density.

           Loops on particles, not on cells.

           :param species_name: The string name of a kinetic-particle species
           :returns: None

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        psaCI = self.pseg_arr[species_name] # segmented array for this species
        (npSeg, pseg) = psaCI.init_out_loop()

        while isinstance(pseg, np_m.ndarray):
            for p in pseg:
                self.number_density_dict[species_name].interpolate_delta_function_to_dofs(p)
            (npSeg, pseg) = psaCI.get_next_segment('out')

        return
#    def accumulate_number_density(self, species_name):ENDDEF

# Useless function:

#class Particle_C(object):
    def set_number_density(self, species_name, value):
        """Set the values in the number density array for the given species.

           :param value: The number density is set to this value.

           :returns: None

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        self.number_density_dict[species_name].set_to_value(value)

        return
#    def set_number_density(self, species_name, value):ENDDEF

#class Particle_C(object):
    def move_particles_in_electrostatic_field(self, dt, neg_electric_field):
        """Apply interpolated electric force to particles.  Compute the
           change in velocity and position in time dt. Use an explicit
           method to integrate the orbit.

           :param str dt: time interval.
           :param neg_electric_field: A Field_C object containing the vector -E.

           :cvar pDim: Number of spatial coordinates in the particle location.

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

#        meshCI = neg_electric_field.meshCI
        pmesh_M = self.pmesh_M

        # Scratch space
        pCoord2 = self.pcoord2 # x,y,z, x0,y0,z0 (or subset)
        dx = self.dx # dx is the distance moved in one step.
#        p_arr = self.one_particle_arr[0]

        pDim = self.particle_dimension

        for sn in self.explicit_species:

            # Invariant parameters
#            print 'sn = ', sn
            qmdt = self.qom[sn]*dt

#        nlocal = self.particles.get_species_particle_count(sp, echoFlag = False)
            psaCI = self.pseg_arr[sn] # segmented array for this species

            if self.get_species_particle_count(sn) == 0: continue

            # negE is long enough to hold one particle segment worth
            # of values.  If segment length is the same for all
            # species, this could be a Field_Particle class array,
            # i.e., one that persists.

#           Accelerate and move all the particles in this species
#            (npSeg, psegIn, psegOut) = psaCI.init_inout_loop()
#            (npSeg, psegIn) = psaCI.init_inout_loop()
            (npSeg, psegIn, psegOut) = psaCI.init_inout_loop()
            particleCount = 0
            # ipOut counts through the 'out' array
            ipOut = 0
            indxChange = False # Flag to indicate if particle SA indices have
                               # changed.  This affects, e.g., trajectories,
                               # which use SA indices to identify particles
                               # chosen for trajectory plots.
            while isinstance(psegIn, np_m.ndarray):
#            print pseg['z']
#                print 'particles_mod: particles module: pseg = ', pseg
                # Compute electric field for each particle

#                print "Before:"
#                print "position:", pseg['x'], pseg['y'] #, pseg['z']
#                print "velocity:", pseg['ux'], pseg['uy'] #, pseg['uz']

                # Interpolate the electric field to the particle positions
                # This is the negative electric field
                neg_electric_field.interpolate_field_to_points(psegIn, self.negE)

                # Truncate negE to the number of particles to get the
                # += operations below to work: These operations need
                # the vectors to be the same length.
#                Eseg = self.negE[0:psegIn.size]
                Eseg = self.negE[0:npSeg]
                # Flip the sign to get correct E. (This is a vector
                # operation on each component of E.)
                for n in Eseg.dtype.names:
                    Eseg[n] *= -1.0


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

                for ipIn in xrange(npSeg):
#                for ip in xrange(pseg.size):

                    # pseg[i] is 'x', 'y', 'z', 'x0',..., 'ux', 'uy',... values of ith item
                    # So pseg[i][0:3] is 'x', 'y', 'z'.
                    # Can't use slice syntax here, because the array data are not homogeneous.

                    # skip deleted particles
                    if psegIn[ipIn]['bitflags'] & self.DELETE_FLAG != 0:
                        indxChange = True # Particle SA indices are stale past
                                          # this point due to deletions.
                        continue

                    # Get the next 'out' segment if the current 'out' segment is
                    # full.
                    if ipOut == self.SEGMENT_LENGTH:
                        psegOut = psaCI.get_next_out_segment()
                        ipOut = 0 # Reset the slot counter for the new segment

#                    if ipOut == 0: psegOut = psaCI.get_next_out_segment()

                    # Copy the 'in' particle to a tmp
#                    p_arr = psegIn[ipIn]
#                    print 'p_arr = ', p_arr
                    psegOut[ipOut] = psegIn[ipIn] # This copies everything over, to ensure weights, flags, etc. get copied.
#                    print 'psegOut = ', psegOut[ipOut]
#                    print 'same?', psegOut[ipOut] is psegIn[ipIn]

                    ## If particle indexing is disturbed, update things that depend on it.
                    # E.g., if this is a trajectory particle, update its SA index
                    if self.traj_T is not None:
                        if psegOut[ipOut]['bitflags'] & self.TRAJECTORY_FLAG != 0:
                            if indxChange is True:
                            # Update the index in the trajectectory
                                self.update_trajectory_particleId(sn, ipIn, ipOut)

                    # Accelerate the particle with the field components in Eseg
                    for comp in Eseg.dtype.names:
                        ucomp = 'u'+comp
#                        psegOut[ipOut][ucomp] = psegIn[ipIn][ucomp] + qmdt*Eseg[ipIn][comp]
                        psegOut[ipOut][ucomp] += qmdt*Eseg[ipIn][comp] # Index OK?

# Instead of this, could do ifs on the dimension of Eseg
                    """
                    psegOut[ipOut]['ux'] = psegIn[ipIn]['ux'] + qmdt*Eseg[ipIn]['x']
                    psegOut[ipOut]['uy'] = psegIn[ipIn]['uy'] + qmdt*Eseg[ipIn]['y']
                    psegOut[ipOut]['uz'] = psegIn[ipIn]['uz'] + qmdt*Eseg[ipIn]['z']
                    """

# Need the dim value

                    # Move the particle
                    # NON-RELATIVISTIC: v and u are the same

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

                    # Check if the particle is still on the meshed region

#                    print 'ip, index =', ip, pseg[ip]['cell_index']
                    pCellIndex = psegOut[ipOut]['cell_index']

#                    print fncName, ": ip, pindex", ip, pCellIndex, "cell index:", pmesh_M.compute_cell_index(pseg[ip])
                    mLastFacet = pmesh_M.NO_FACET
                    facetCrossCount = 0
                    while not pmesh_M.is_inside(psegOut[ipOut], pCellIndex):
                        # The particle has left this cell.  We
                        # need to track it across each facet in case
                        # there's a boundary-condition on that facet.
#                        print fncName, "particle has migrated"

                        facetCrossCount += 1
                        # Check for an abnormal number of facet crossings:
                        if facetCrossCount > self.MAX_FACET_CROSS_COUNT:
                            errorMsg = "%s !!! MAX_FACET_CROSS_COUNT exceeded!!!" % (fncName)
                            sys.exit(errorMsg)

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

                        dx = pCoord2[0:pDim] - pCoord2[pDim:2*pDim] # Move vector

# See save10/Particle_Module.py for a different search, where only nearby cells are looked at.
                        (cFacet, pathFraction) = pmesh_M.find_facet(pCoord2[pDim:2*pDim], dx, pCellIndex)
                        if cFacet != pmesh_M.NO_FACET:

#                        if found_cell != True: print "Particle is not in nearby cells; using BB search"
#                        ci = pmesh_M.compute_cell_index(pseg[ip])
#                        print "Found particle is in cell", ci
#                        pseg[ip]['cell_index'] = ci

                            # Compute the crossing point
#                            print fncName, "Before truncation p =", psegOut[ipOut]
                            i=0
                            for coord in self.position_coordinates:
                                coord0 = coord+'0'
                                psegOut[ipOut][coord0] = psegOut[ipOut][coord0] + pathFraction*dx[i]
                                i+=1
#                            print "After truncation p =", psegOut[ipOut]

                            # Look up the mesh-level index of this facet...
                            mFacet = pmesh_M.cell_entity_index_dict['facet'][pCellIndex][cFacet]
# mFacet should never be the same as the last facet crossed: check this
                            if mFacet == mLastFacet:
                                errorMsg = "%s The mesh index of the facet crossed is %d, the same as the last facet crossed. This should not happen since the particle cannot cross the same facet twice in one move!" % (fncName, mFacet)
                                sys.exit(errorMsg)
                            else:
                                mLastFacet = mFacet
                            # ...and get the value of the facet marker.
# mFacet is a numpy.uint32 (size_t), but the FacetFunction wants an int argument.
#                            print "type is:", type(mFacet)
                            facValue = pmesh_M.particle_boundary_marker[mFacet]
                            # Check if this facet has a non-zero marker
                            if facValue != 0:
                                # Call the function associated with this value.

# Maybe this needs the whole particle object?  It could be the end-of-the line for the particle.
# Need to check what needs to be in the following arg list:
                                self.pmesh_bcCI.bc_function_dict[facValue][sn](psegOut[ipOut], sn, mFacet)

# Not sure if this else is needed:                            else:

                            # Update the cell index to the new cell
                            pCellIndex = pmesh_M.cell_neighbor_dict[pCellIndex][cFacet]
                            psegOut[ipOut]['cell_index'] = pCellIndex
                            # If the particle has left the grid, exit this 'while' loop.
                            if pCellIndex == pmesh_M.NO_CELL: break # Exit the 'while' loop
                            # If pathFraction is zero, this is a flag
                            # to verify that the index is correct
#                            if pathFraction == 0.0:
                        else:
                            errorMsg = "%s The cell index of the facet crossed is %d. This should not happen since the particle has left its initial cell cell!" % (fncName, cFacet)
                            sys.exit(errorMsg)
#                       END:if cFacet != pmesh_M.NO_FACET:
#
#                   END:while not pmesh_M.is_inside(psegOut[ipOut], pCellIndex)

                    # Don't need this since we just look up the cell
                    # index when computing negE above
                    # pseg[ip]['cell_index'] = compute_cell_index(Point(p))

                    # Check that this particle has not been deleted before
                    # incrementing the 'out' particle counter
                    if psegOut[ipOut]['bitflags'] & self.DELETE_FLAG == 0:
                        ipOut += 1
                    else:
#                        print "Particle", ipIn, "of species", sn, "has been deleted."
                        # Remove the particle from the trajectory list if it was tagged
                        if self.traj_T is not None:
                            if psegIn[ipIn]['bitflags'] & self.TRAJECTORY_FLAG != 0:
                                self.remove_trajectory_particleId(sn, ipIn, psegOut[ipOut])

                    # Check if we've reached the end of this segment.  If
                    # so, we need to start writing on a new segment.
                    if (ipOut == self.SEGMENT_LENGTH):
                        particleCount += self.SEGMENT_LENGTH
#                        ipOut = 0 # This will cause get_next_out_segment() to be called
                                   # above if there are more 'in' particles to be processed.

                # Done with this segment.
                # Get the next one, if it exists.
                (npSeg, psegIn) = psaCI.get_next_segment('in')
#                pseg = psaCI.get_next_segment()
                # Loop over particles in this segment ends

            # Set values that will be used to keep track of the particle arrays
            particleCount += ipOut
            psaCI.set_number_of_items('out', particleCount)
            # Loop over segmented array ends

# Compute new density here?

        # Loop over Species ends
        return
#    def move_particles_in_electrostatic_field(self, dt, neg_electric_field):ENDDEF

#class Particle_C(object):
    def move_neutral_particles(self, dt):
        """Move neutral particles on a mesh.

           Compute change in position in time dt. Use an explicit method to integrate
           the orbit.

           Arguments:
               dt: time interval.

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        # Scratch space
        pCoord2 = self.pcoord2 # x,y,z, x0,y0,z0 (or subset)
        dx = self.dx # dx is the distance moved in one step.
#        p_arr = self.one_particle_arr[0]

        pmesh_M = self.pmesh_M
        pDim = self.particle_dimension

        for sn in self.neutral_species:

            # Invariant parameters
#            print 'sn = ', sn

#        nlocal = self.particles.get_species_particle_count(sp, echoFlag = False)
            psaCI = self.pseg_arr[sn] # segmented array for this species

            if self.get_species_particle_count(sn) == 0: continue

#           Move all the particles in this species
            (npSeg, psegIn, psegOut) = psaCI.init_inout_loop() # (number of particles in
#            (npSeg, psegIn) = psaCI.init_inout_loop() # (number of particles in
                                                      # this segment, ref to
                                                      # segment)
            particleCount = 0
            # ipOut counts particles being written to the current
            # 'out' segment.
            ipOut = 0
            indxChange = False # Flag to indicate if particle SA indices have
                               # changed.  This affects, e.g., trajectories,
                               # which use SA indices to identify particles
                               # chosen for trajectory plots.
            while isinstance(psegIn, np_m.ndarray): # Keep looping until we run
                                                    # out of 'in' segments
                for ipIn in xrange(npSeg): # Loop on the particles in this 'in'
                                           # segment
                    # psegIn[ipIn] has the 'x', 'y', 'z', 'x0',..., 'ux', 'uy',... values of ith item
                    # So psegIn[ipIn][0:3] is 'x', 'y', 'z'.
                    # Can't use slice syntax here, because the array data are not of homogeneous type.

                    # Skip deleted particles
                    if psegIn[ipIn]['bitflags'] & self.DELETE_FLAG != 0:
                        indxChange = True # Particle SA indices are stale past
                                          # this point due to deletions.
                        continue # skip to the next particle

                    # If this 'out segment is full, get the next 'out' segment,
                    # provided there are still some 'in' particles to advance.
                    if ipOut == self.SEGMENT_LENGTH:
                        psegOut = psaCI.get_next_out_segment()
                        ipOut = 0 # Reset the counter for the new segment
#                    if ipOut == 0: psegOut = psaCI.get_next_out_segment()

                    # The following COPIES data, rather than just setting a
                    # reference to existing data.
                    psegOut[ipOut] = psegIn[ipIn] # Copy this particle's data
                                                  # from the input slot to the
                                                  # output slot

# Also, check if the following copies, or resets the reference. Ans: it COPIES.
#                    p_out = psegOut[ipOut]
# To get a REFERENCE, need to use SLICING syntax:
#                    p_out = psegIn[ipIn][:] # This copies everything over, to ensure weights, flags, etc. get copied.
# But you can't get a reference to a single element of psegIn[]. See:
# http://stackoverflow.com/questions/23654088/reference-of-a-single-numpy-array-element
# See also HPL p. 137. b=a[1,:], but b is still an array.
# What about getting a ref to one particle returned from a function? That's possible because it uses the stack?

                    ## If particle indexing is disturbed, update things that depend on it.

                    # E.g., if this is a trajectory particle, update its SA index
                    if self.traj_T is not None:
                        if psegOut[ipOut]['bitflags'] & self.TRAJECTORY_FLAG != 0:
                            if indxChange is True:
                            # Update the index in the trajectectory
                                self.update_trajectory_particleId(sn, ipIn, ipOut)

                    ## Move the particle
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
                    mLastFacet = pmesh_M.NO_FACET
                    facetCrossCount = 0
                    while not pmesh_M.is_inside(psegOut[ipOut], pCellIndex):
                        # The particle has left this cell.  We
                        # need to track it across each facet in case
                        # there's a boundary-condition on that facet.
#                        print fncName, "particle has migrated"

                        facetCrossCount += 1
                        # Check for an abnormal number of facet crossings:
                        if facetCrossCount > self.MAX_FACET_CROSS_COUNT:
                            errorMsg = "%s !!! MAX_FACET_CROSS_COUNT exceeded!!!" % (fncName)
                            sys.exit(errorMsg)

                        # Compute dx[], the move vector that starts in
                        # the current cell
                        i=0
                        for coord in self.position_coordinates:
                            pCoord2[i] = psegOut[ipOut][coord] # Present position
                            coord0 = coord+'0'
                            pCoord2[i+pDim] = psegOut[ipOut][coord0] # Start of path in this cell
                            i+=1
#                        pCoord2[:] = psegOut[ipOut] # Alias to the position coordinates
                        dx = pCoord2[0:pDim] - pCoord2[pDim:2*pDim] # Move vector

# could return the crossing-point in pCoord2[]
                        (cFacet, pathFraction) = pmesh_M.find_facet(pCoord2[pDim:2*pDim], dx, pCellIndex)
                        if cFacet != pmesh_M.NO_FACET:
#                            print "facet crossed is", cFacet

                            # Compute the crossing point
#                           print fncName, "Before truncation p =", psegOut[ipOut]
                            i=0
                            for coord in self.position_coordinates:
                                coord0 = coord+'0'
                                psegOut[ipOut][coord0] = psegOut[ipOut][coord0] + pathFraction*dx[i]
                                i+=1
#                            print "After truncation p =", psegOut[ipOut]

                            # Look up the mesh-level index of this facet...
                            mFacet = pmesh_M.cell_entity_index_dict['facet'][pCellIndex][cFacet]
# mFacet should never be the same as the last facet crossed: check this
                            if mFacet == mLastFacet:
                                errorMsg = "%s The mesh index of the facet crossed is %d, the same as the last facet crossed. This should not happen since the particle cannot cross the same facet twice in one move!" % (fncName, mFacet)
                                sys.exit(errorMsg)
                            else:
                                mLastFacet = mFacet
                            # ...and get the value of the facet marker.
# mFacet is a of type 'size_t' or numpy.uint32, but the FacetFunction wants an int argument.
#                            print "type is:", type(mFacet)
                            facValue = pmesh_M.particle_boundary_marker[mFacet]
                            # Check if this facet has a non-zero marker
                            if facValue != 0:
                                # Call the function associated with this value.

# Maybe this needs the whole particle object?  It could be the end-of-the line for the particle.
# Need to check what needs to be in the following arg list:
                                self.pmesh_bcCI.bc_function_dict[facValue][sn](psegOut[ipOut], sn, mFacet)

# Not sure if this else is needed:                            else:

#                            # Update the cell index to the new cell
#                                pCellIndex = pmesh_M.cell_neighbor_dict[pCellIndex][facet]
#                                psegOut[ipOut]['cell_index'] = pCellIndex
#                                print "cell_index updated to", psegOut[ipOut]['cell_index']
                            # Update the cell index to the new cell
                            pCellIndex = pmesh_M.cell_neighbor_dict[pCellIndex][cFacet]
                            psegOut[ipOut]['cell_index'] = pCellIndex
                            # If the particle has left the grid, exit this 'while' loop.
                            if pCellIndex == pmesh_M.NO_CELL: break # Exit the 'while' loop

                        else:
                            errorMsg = "%s The cell index of the facet crossed is %d. This should not happen since the particle has left its initial cell cell!" % (fncName, cFacet)
                            sys.exit(errorMsg)
#                       END:if cFacet != pmesh_M.NO_FACET:
#                   END:while not pmesh_M.is_inside(psegOut[ipOut], pCellIndex)

                    # Don't need this since we just look up the cell
                    # index when computing negE above
                    # pseg[ip]['cell_index'] = compute_cell_index(Point(p))

                    # Check that this particle has not been deleted before
                    # incrementing the counter
                    if psegOut[ipOut]['bitflags'] & self.DELETE_FLAG == 0:
                        ipOut += 1
                    else:
#                        print "Particle", ipIn, "of species", sn, "has been deleted."
                        pass

                    # Check if we've reached the end of this segment.  If
                    # so, we need to start writing on a new segment.
                    if (ipOut == self.SEGMENT_LENGTH):
                        particleCount += self.SEGMENT_LENGTH
#                        ipOut = 0 # This will cause get_next_out_segment() to be called
                                   # above, if there are more 'in' particles to be processed.

                # Done with this segment.
                # Get the next one, if it exists.
                (npSeg, psegIn) = psaCI.get_next_segment('in')
#                pseg = psaCI.get_next_segment()
                # Loop over particles in this segment ends

            # Set values that will be used to keep track of the particle arrays
            particleCount += ipOut
            psaCI.set_number_of_items('out', particleCount)
            # Loop over segmented array ends

# Compute new density here?

        # Loop over Species ends
        return
#    def move_neutral_particles(self, dt, neg_electric_field):ENDDEF

#class Particle_C(object):
    def move_particles_in_electrostatic_potential(self, species_names, dt, fCI, fpCI):
        """Do an energy-conserving push in the electrostatic
           potential.
        """
        for sn in species_names:

            # Invariant parameters
            qmdt = self.qom[sn]*dt

#        nlocal = self.particles.get_species_particle_count(sp, echoFlag = False)
            psaCI = self.pseg_arr[sn] # segmented array for this species

            # Eseg is long enough to hold one particle segment worth
            # of values.  If segment length is the same for all
            # species, this could be a Field_Particle class array,
            # i.e., one that persists.

#            Eseg = np.empty(len(psaCI[0]), dtype=fpCI.Eseg_dict)

#           Accelerate and move all the particles in this species
            psaCI.init_segment_loop()
            pseg = psaCI.get_next_segment()
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
                pseg = psaCI.get_next_segment()

        return
#    def move_particles_in_electrostatic_potential(self, species_names, dt, fCI, fpCI):ENDDEF



#class Particle_C(object):
    def move_particles_in_uniform_fields(self, species_name, ctrl, print_flag = False):
        """Apply electric field ctrl.E0 to particles.  Compute the
           resulting change in particle velocities and positions in time
           ctrl.dt.
        """
# Use numpy array syntax to do this, instead of loops; see Susie; See HPL Sec 4.2.2

        dt = ctrl.dt
        E0 = ctrl.E0

        # Invariant parameters

        qmdt = self.qom[species_name]*dt

#        nlocal = self.particles.get_species_particle_count(sp, echoFlag = False)

        # Accelerate the particles
        psaCI = self.pseg_arr[species_name] # segmented array for this species

        # This loop moves particles, so swap the in/out particle arrays
#        (npSeg, psegIn, psegOut) = psaCI.init_inout_loop()
#        (npSeg, psegIn) = psaCI.init_inout_loop()
        (npSeg, psegIn, psegOut) = psaCI.init_inout_loop()

#        print "move_particles_in_uniform_fields: npSeg, psegIn, psegOut:", npSeg, psegIn, psegOut
        # Get the first segment of particles to move
#        pseg = psaCI.get_next_segment()

        # Loop on segments until a None value is returned
        particleCount = 0
        # ipOut counts through the 'out' array
        ipOut = 0
        while isinstance(psegIn, np_m.ndarray):
#            print "type pseg = ", type(pseg)
#            print pseg['z']
#            print pseg
            # Accelerate this block of particles

            # Loop on the particles in this segment.
            # If a particle has been deleted, skip it.

            # ipIn counts through the 'in' segment
            for ipIn in xrange(npSeg):

                ## Skip deleted particles
                if psegIn[ipIn]['bitflags'] & self.DELETE_FLAG != 0:
#                    print "Deleting particle with ipIn =", ipIn
                    indxChange = True # Particle SA indices are stale past
                                      # this point due to deleted items.
                    continue # skip to the next particle

                # Only gets the next 'out' segment if there are more 'in' particles to do

                if ipOut == self.SEGMENT_LENGTH:
                    psegOut = psaCI.get_next_out_segment()
                    ipOut = 0 # Reset the counter for the new segment
#                if ipOut == 0: psegOut = psaCI.get_next_out_segment()

                psegOut[ipOut] = psegIn[ipIn] # This copies everything over, to ensure weights, flags, etc. get copied.

                ## If particle indexing is disturbed, update things that depend on it.

                # E.g., if this is a trajectory particle, update its SA index
                if self.traj_T is not None:
                    if psegOut[ipOut]['bitflags'] & self.TRAJECTORY_FLAG != 0:
                        if indxChange is True:
                        # Update the index in the trajectectory
                            (full_index_in, full_index_out) = self.update_trajectory_particleId(species_name, ipIn, ipOut)
#                            print "Updated trajectory index", ipIn, "to", ipOut
                            print "Updated trajectory index", full_index_in, "to", full_index_out

                # Accelerate the particle
                psegOut[ipOut]['ux'] = psegIn[ipIn]['ux'] + qmdt*E0.x
                psegOut[ipOut]['uy'] = psegIn[ipIn]['uy'] + qmdt*E0.y
                psegOut[ipOut]['uz'] = psegIn[ipIn]['uz'] + qmdt*E0.z

                # Move the particle
                # NON-RELATIVISTIC: v and u are the same
                psegOut[ipOut]['x'] = psegIn[ipIn]['x'] + psegOut[ipOut]['ux']*dt
                psegOut[ipOut]['y'] = psegIn[ipIn]['y'] + psegOut[ipOut]['uy']*dt
                psegOut[ipOut]['z'] = psegIn[ipIn]['z'] + psegOut[ipOut]['uz']*dt

                ipOut += 1
                # Check if we've reached the end of this segment.  If
                # so, we need to start writing on a new segment.  But
                # don't increment to a new 'out' segment unless there
                # are actually more 'in' particles to process.
                if (ipOut == self.SEGMENT_LENGTH):
                    particleCount += ipOut

    #            print pseg['z']
            (npSeg, psegIn) = psaCI.get_next_segment('in')

        # Set values that will be used to keep track of the particle arrays
#        print "Species", species_name, "npSeg=", npSeg, "psegIn", psegIn, "particleCount", particleCount
        particleCount += ipOut
        psaCI.set_number_of_items('out', particleCount)

        return
#    def move_particles_in_uniform_fields(self, species_name, ctrl, print_flag = False):ENDDEF

# Just push one segment:

    def move_particle_segment_in_uniform_fields(self, dt, qmdt, pseg, E0):
        """Push one segment of particles in fields that are uniform in space and
           time.
        """
# Use numpy array syntax to do this, instead of loops; see Susie; See HPL Sec 4.2.2

            # Accelerate this block of particles
        pseg['ux'] += qmdt*E0.x
        pseg['uy'] += qmdt*E0.y
        pseg['uz'] += qmdt*E0.z
# NON-RELATIVISTIC: v and u are the same
            # Move the particles
        pseg['x'] += pseg['ux']*dt
        pseg['y'] += pseg['uy']*dt
        pseg['z'] += pseg['uz']*dt
#            print pseg['z']

        return
#    def move_particles_in_uniform_fields(self, species_name, ctrl, print_flag = False):ENDDEF

#class Particle_C(object):
    def update_trajectory_particleId(self, sn, i_in, i_out):
        """Replace the full index of a trajectory particle with the new value.

           Trajectory particles are identified by their full SA index.  When
           particles are deleted from the SA, these indices change for the
           remaining particles, and have to be updated.

           :param str sn: The name of the current species being advanced.
           :param int i_in: The particle's index in the current 'in' segment.
           :param int i_out: The particle's index in the current 'out' segment.
        """

        psaCI = self.pseg_arr[sn] # The segmented array for this species

        # Obtain the full indices of this particle in the 'in' and 'out' arrays.
        (full_index_in, full_index_out) = psaCI.get_full_indices(i_in, i_out)

        # Find the position of full_index_in:
        indx = self.traj_T.ParticleIdList[sn].index(full_index_in)
        # Replace this with the new full index
        self.traj_T.ParticleIdList[sn][indx] = full_index_out

        return full_index_in, full_index_out
#    def update_trajectory_particleId(self, sn, i_in, i_out):ENDDEF

#class Particle_C(object):
    def remove_trajectory_particleId(self, sn, i_in, p, step, time):
        """Record the last data for a particle and remove it from the trajectory
           list.

           When a particle that is tagged as a trajectory particle is deleted,
           it has to be removed for the trajectory list.  Record the last datum
           for this particle, and then change the particle index to NO_PINDEX to
           skip it in subsequent looping over trajectory particles.

           :param str sn: The name of the current species being advanced.
           :param int i_in: The particle's index in the current 'in' segment.
           :param p: The current attributes of the particle being deleted.

           :cvar NO_PINDEX: a flag indicating the particle no longer exists.
        """

        traj_T = self.traj_T
        psaCI = self.pseg_arr[sn] # The segmented array for this species
        p_arr = self.one_particle_arr

        # Obtain the full index of this particle
        full_index_in = psaCI.get_full_index(i_in, 'in')

        # Record the last position
        # Retrieve particle using its full index
#        p_arr[0] = psaCI.get(full_index_in)
        p_arr[0] = p

        # Copy the particle values into the trajectory

HERE

        # Find the position of full_index_in:
        indx = self.traj_T.ParticleIdList[sn].index(full_index_in)

        count = self.traj_T.TrajectoryLength[sn][indx]
        for comp in traj_T.explicit_dict['names']:
            if comp in p_arr.dtype.names:
                traj_T.DataList[sn][indx][comp][count] = p_arr[0][comp]
                # The particle has been deleted and may be out-of-bounds, so set
                # E to zero instead of interpolating.
                # This isn't right: 'Ex' never passes this test.
            if comp in self.negE.dtype.names:
                Ecomp = 'E'+comp
                traj_T.DataList[sn][indx][Ecomp][count] = 0.0

        self.traj_T.TrajectoryLength[sn][indx] += 1 # Increment the trajectory length

        # Mark this as a trajectory where the particle no longer exists.
        self.traj_T.ParticleIdList[sn][indx] = self.traj_T.NO_PINDEX

        return full_index_in
#    def remove_trajectory_particleId(self, sn, i_in):ENDDEF

#class Particle_C(object):
    def record_trajectory_datum(self, species_name, p, step, time, neg_E_field):
        """Save a single data-record of a particle trajectory.

           :param species_name: Name of the species that the marked particle belongs
                                to.
           :param p: One particle record, as defined by particle_dtype in
                     Particle_C.__init__
           :param neg_E_field: A Field_C object containing the vector -E.

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        traj_T = self.traj_T
        p_arr = self.one_particle_arr

        # Copy the particle into an array, since an array must be passed below.
        p_arr[0] = p

        # print "record_trajectory_datum: p_arr[0] = ", p_arr[0]

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

        # print "E_arr.dtype.names =", E_arr.dtype.names

        # Copy the particle values into the new trajectory, which is the last one
        # added onto the list of trajectories for this species.

        inew = len(traj_T.ParticleIdList[species_name]) - 1 #
        count = traj_T.TrajectoryLength[species_name][inew]
        for comp in traj_T.explicit_dict['names']:
            # print 'comp = ', comp
            if comp == 'step':
                traj_T.DataList[species_name][inew][comp][count] = step
            elif comp == 't':
                traj_T.DataList[species_name][inew][comp][count] = time
            elif comp in p_arr.dtype.names:
                traj_T.DataList[species_name][inew][comp][count] = p_arr[0][comp]
            elif comp in E_arr.dtype.names:
                traj_T.DataList[species_name][inew][comp][count] = E_arr[0][comp]

        traj_T.TrajectoryLength[species_name][inew] += 1 # Increment the trajectory length

        return
#    def record_trajectory_datum(self, neg_E_field=None):ENDDEF

#class Particle_C(object):
    def record_trajectory_data(self, step, time, neg_E_field=None):
        """Save one trajectory data-record for all the marked particles of all
           species.

           :param int step: The current timestep
           :param float time: The current physical time
           :param neg_E_field: A Field_C object containing the vector -E

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        # Avoid a second call if the timestep hasn't changed
        if step == self.traj_T.last_step:
            print fncName, "(DnT INFO) Second call on timestep %d: return without storing the data." % step
            return

        traj_T = self.traj_T
        p_arr = self.one_particle_arr
#        print 'p_arr.dtype.names = ', p_arr.dtype.names
#        count = self.traj_T.count

        print 'step = ', step

        # Loop on trajectories of explicit particles
        for sp in traj_T.explicit_species:
            pseg_arrCI = self.pseg_arr[sp] # The SA for this species
            for i in xrange(len(traj_T.ParticleIdList[sp])): # i loops on
                                                             # particle-index array
#            for ip in traj_T.ParticleIdList[sp]:
                ip = traj_T.ParticleIdList[sp][i] # Look up the full particle index

                # If a particle no longer exists, skip
                if ip == traj_T.NO_PINDEX: continue

                # Retrieve particle using its full index
                p_arr[0] = pseg_arrCI.get(ip) # pulls value from the 'out' array.

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

                # Copy the particle values into the trajectory
                count = traj_T.TrajectoryLength[sp][i]
                for comp in traj_T.explicit_dict['names']:
#                    print 'comp = ', comp
                    if comp == 'step':
                        traj_T.DataList[sp][i][comp][count] = step
                    elif comp == 't':
                        traj_T.DataList[sp][i][comp][count] = time
                    elif comp in p_arr.dtype.names:
                        traj_T.DataList[sp][i][comp][count] = p_arr[0][comp]
                    elif comp in E_arr.dtype.names:
                        traj_T.DataList[sp][i][comp][count] = E_arr[0][comp]

#        print 'traj_T.DataList: ', traj_T.DataList['testelectrons'][0]

                traj_T.TrajectoryLength[sp][i] += 1 # Increment the trajectory length

        # Loop on trajectories of implicit particles

        # Update "last_step" before returning
        self.traj_T.last_step = step

        return
#    def record_trajectory_data(self, neg_E_field=None):ENDDEF


#class Particle_C(object):
    def accumulate_charge_density_from_particles(self, charge_density):
        """Compute the charge density contributed by the kinetic particles.

           :param charge_density: Storage into which charge-density is accumulated.

           :returns: None

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        ## Loop over species
        if self.pmesh_M is not None:
            for s in self.species_names:
                charge = self.charge[s]
                # Loop over the number-density arrays
                charge_density.multiply_add(self.number_density_dict[s], charge)
        else:
            print fncName, "(DnT WARNING) The reference to pmesh_M is None"

        return
#    def accumulate_charge_density_from_particles(self, charge_density):ENDDEF


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
        pseg_arrCI = self.pseg_arr[speciesName] # The SegmentedArray_C object for this species
        pDim = self.particle_dimension

        # References to scratch space
        particle = self.one_particle_arr[0]
        pCoord = self.pcoord # x, y, z (or subset)
        pVel = self.pvel # ux, uy, uz (or subset)
        random_vals = self.random_vals

        # Loop over the domain creating particles in each cell

        pCoord[:] = 0.0 # Zero out this array for the case where there are more
                        # particle spatial dimensions than mesh spatial
                        # dimensions

        # Multiplier for particle weight
        weightMult = numberDensity/numberPerCell
        if self.coordinate_system == '1D-spherical-radius':
            weightMult /= 4.0*np_m.pi # 4 \pi radians in a sphere

        for icell in xrange(nCell):
            # Compute the particle weight (number of particles per macroparticle)
            weight = domain.volume[icell]*weightMult

            cellIndex = domain.cell_index[icell]

            # Set positions and velocities
            cellRad = domain.radius[icell]
            cellMid = domain.midpoint[icell]
            gDim = domain.gdim

            for ip in xrange(numberPerCell):
                # Put the first particle at the centroid
                # Put the rest at random positions in the cell
                bitflags = 0b00 # bit flags variable is all zeroes
                if ip == 0:
                    for i in range(gDim):
                        pCoord[i] = cellMid[i]
                        # Turn trajectory flag ON for the first particle:
                        bitflags = bitflags | Particle_C.TRAJECTORY_FLAG
                else:
                    while True:
                        # Put the particle in a sphere
                        random_vals = np_m.random.uniform(-1.0, 1.0, gDim)
                        for i in range(gDim):
                            pCoord[i] = cellMid[i]+random_vals[i]*cellRad
                        # Check if its in the cell
                        if domain.in_cell(pCoord, icell):
                            break

                # Generate a velocity vector

# For 'r_theta_z' this needs to change?

                pVel = np_m.random.normal(0.0, thermalSpeed, pDim) + driftVelocity

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

                # Store the particle
                p, pindex = pseg_arrCI.put(particle)

                # If the particle is tagged as a trajectory particle, initialize its
                # trajectory information. Note: if this is being called by
                # initialize_particles(), the E-field has not been computed yet.
                if p['bitflags'] & self.TRAJECTORY_FLAG != 0:
                    if self.traj_T is not None:
#                    print 'pindex for trajectory = ', pindex
                        dynamicsType = self.dynamics[speciesName]
                        self.traj_T.create_trajectory(speciesName, pindex, dynamicsType)
                        # Record the initial datum for this new trajectory
                        self.record_trajectory_datum(speciesName, p, step, time, neg_E_field)
                    else:
# Instead of printing this message, a traj_T object could be created here?
                        print fncName, "(DnT WARNING) A trajectory flag is on, but no trajectory object has been created yet."

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
                sys.exit(errorMsg)

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
        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        groupName = "Step#" + str(self.h5_step_counter)
        group = self.particle_output_handle.create_group(groupName)

        ### Count the total number of particles and see if h5_buffer_length needs to
        ### be increased
        totalParticleCount = self.get_total_particle_count(print_flag = False)
#        print fncName, "totalParticleCount =", totalParticleCount, "h5_buffer_length =", self.h5_buffer_length
        if totalParticleCount > self.h5_buffer_length:
            newNumSegs = 1 + totalParticleCount/self.SEGMENT_LENGTH # integer arithmetic
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
                if np_m.issubdtype(self.particle_dtype[pA], np_m.float):
                    h5Buf.dtype = np_m.float64
                elif np_m.issubdtype(self.particle_dtype[pA], np_m.integer):
                    h5Buf.dtype = np_m.int64
                else:
                    errorMsg = fncName + "The type of particle attribute " + pA + " is not float or integer. It is " + str(self.particle_dtype[pA])
                    sys.exit(errorMsg)

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
                    ## Loop over the segments of the 'out' array to get the attribute pA
                    psaCI = self.pseg_arr[s] # segmented array for this species
                    (npSeg, pseg) = psaCI.init_out_loop()
                    particleCount = 0
                    while isinstance(pseg, np_m.ndarray):
#                        print "Array for", pA, "is: ", pseg[pA], "shape = ", pseg.shape
#                        print "Range of h5Buf is:", aOff, aOff+npSeg, "shape =", h5Buf.shape
                        h5Buf[aOff:aOff+npSeg] = pseg[pA]
                        aOff += npSeg
                        (npSeg, pseg) = psaCI.get_next_segment('out')

#            print "h5Buf is:", h5Buf[0:aOff]
            dset = group.create_dataset(pA, data=h5Buf[0:aOff])
#                    group.create_dataset(pA, data=h5Buf[:])
#                    d2 = group.create_dataset("y",data=y[:], dtype='f')

        self.h5_step_counter += 1

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
    def __init__(self, speciesNames, pmesh_M, userParticleMeshFunctionsClass, print_flag = False):
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

        # We want to associate the names of the user-supplied boundary
        # functions with the int values of the facet tags.  First, we
        # swap the BC dictionary keys and values: e.g., 'xmin': 1 to
        # 1: 'xmin'. This lets us use the facet uint (size_t) as an index to
        # get the name of the boundary.

        pBDictInv = {v: k for k, v in pBDict.iteritems()}

        particleBoundaryTags = pBDictInv.keys()

        # Initialize a new dictionary that will contain the name of a boundary
        # function for each [boundary tag][species] pair (i.e., the key is a 2D object).
        self.bc_function_dict = dict((intTag, dict((sp, None) for sp in speciesNames)) for intTag in particleBoundaryTags)

        # Find the global default boundary function, if there is one,
        # in the supplied UserParticleBoundaryFunctions object.
        bcFunctionName = 'default_bc'
        if hasattr(userParticleMeshFunctionsClass, bcFunctionName):
            bcGlobalDefaultFunction = getattr(userParticleMeshFunctionsClass, bcFunctionName)
        else:
            bcGlobalDefaultFunction = None

        # Loop on particle boundary tags and on particle species
        # names to find the most specific BC.

        # Overwrite the global default function with a function specific to
        # each boundary, if there is one.
        for intTag in particleBoundaryTags:
            bcFunctionName = 'default_bc_at_' + pBDictInv[intTag]
            if hasattr(userParticleMeshFunctionsClass, bcFunctionName):
                bcBoundaryDefaultFunction = getattr(userParticleMeshFunctionsClass, bcFunctionName)
            else:
                bcBoundaryDefaultFunction = bcGlobalDefaultFunction
            for sp in speciesNames:
                # Overwrite the default for this boundary with a
                # function specific to this boundary and species, if
                # there is one.
                bcFunctionName = 'bc_at_' + pBDictInv[intTag] + '_for_' + sp
                if hasattr(userParticleMeshFunctionsClass, bcFunctionName):
                    bcFunction = getattr(userParticleMeshFunctionsClass, bcFunctionName)
                else:
                    bcFunction = bcBoundaryDefaultFunction
                if bcFunction is None:
                    print "ParticleMeshBoundaryConditions_C: No boundary condition specified for", pBDictInv[intTag], "/", sp
                elif print_flag:
                    print "ParticleMeshBoundaryConditions_C: Boundary condition for", pBDictInv[intTag], "/", sp, "is", bcFunction
                self.bc_function_dict[intTag][sp] = bcFunction
        return
#    def __init__(self, particleInputCI, particleCI, print_flag = False):ENDDEF

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
    def __init__(self, speciesNames, pmesh_M, userParticleMeshSourceFunctionsClass, print_flag = False):
        """Initialize particle source functions.

        """
        # Set local names from passed parameters
        pSrcDict = pmesh_M.particle_source_dict

        # We want to associate the names of the user-supplied cell functions with the int
        # values of the cell tags.  First, we swap the source dictionary keys and values:
        # e.g., 'xmin': 1 to 1: 'xmin'. This lets us use the cell-tag as an index to get
        # the name of the source.

        pSrcDictInv = {v: k for k, v in pSrcDict.iteritems()}

        particleSourceTags = pSrcDictInv.keys()

        # Initialize a new dictionary that will contain the name of a
        # boundary function for each (boundary, species) pair.
        # The first index is the integer tag that
        # marks the boundary facets.
        self.srcFunctionDict = dict((intTag, dict((sp, None) for sp in speciesNames)) for intTag in particleSourceTags)

        # Find the global default boundary function, if there is one,
        # in the supplied UserParticleBoundaryFunctions object.
        srcFunctionName = 'default_bc'
        if hasattr(userParticleMeshSourceFunctionsClass, srcFunctionName):
            bcGlobalDefaultFunction = getattr(userParticleMeshSourceFunctionsClass, srcFunctionName)
        else:
            bcGlobalDefaultFunction = None

        # Loop on particle boundary tags and on particle species
        # names to find the most specific BC.

        # Overwrite the global default function with a function specific to
        # each boundary, if there is one.
        for intTag in particleSourceTags:
            srcFunctionName = 'default_bc_at_' + pSrcDictInv[intTag]
            if hasattr(userParticleMeshSourceFunctionsClass, srcFunctionName):
                bcBoundaryDefaultFunction = getattr(userParticleMeshSourceFunctionsClass, srcFunctionName)
            else:
                bcBoundaryDefaultFunction = bcGlobalDefaultFunction
            for sp in speciesNames:
                # Overwrite the default for this boundary with a
                # function specific to this boundary and species, if
                # there is one.
                srcFunctionName = 'bc_at_' + pSrcDictInv[intTag] + '_for_' + sp
                if hasattr(userParticleMeshSourceFunctionsClass, srcFunctionName):
                    srcFunction = getattr(userParticleMeshSourceFunctionsClass, srcFunctionName)
                else:
                    srcFunction = bcBoundaryDefaultFunction
                if srcFunction is None:
                    print "ParticleMeshSources_C: No boundary condition specified for", pSrcDictInv[intTag], "/", sp
                elif printTag:
                    print "ParticleMeshSources_C: Boundary condition for", pSrcDictInv[intTag], "/", sp, "is", srcFunction
                self.srcFunctionDict[intTag][sp] = srcFunction

        return
#    def __init__(self, particleInputCI, particleCI, print_flag = False):ENDDEF

#class ParticleMeshSources_C(object):ENDCLASS
