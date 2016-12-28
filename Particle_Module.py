# Particle module

"""Particle_Module treats discrete particles.
"""

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
__all__ = ['Particle_C.initial_distribution_type', 
           'Particle_C.number_of_species', 
           'Particle_C.index', 
           'Particle_C.name',
           'Particle_C.number_per_cell', ]



import sys
import numpy as np_M

# use initial underscores for locals if you want to use a non-local name

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

       :cvar position_coordinates: Labels of the position coordinate axes.
       :vartype position_coordinates: string array
       :cvar dimension: Number of spatial coordinates in a particle position.

    """

    # Static class variables

    # Particles are stored in segments of the following length:
    SEGMENT_LENGTH = 100
#    SEGMENT_LENGTH = 5
    # Bitmasks are static class variables
    DELETE_FLAG = 0b01 # the lowest bit is 1
    TRAJECTORY_FLAG = 0b10 # the second lowest bit is 1

#    def __init__(self, species_input, phase_coordinates, precision, segment_length, user_particles_class, echoFlag):


#    def __init__(self, species_input, phase_coordinates, precision, user_particles_class, echoFlag):
#class Particle_C(object):
    def __init__(self, particleInputCI, printFlag = False):
        """Take a list of kinetic species provided by the user and create the initial plasma
        """
        # Set local variables from passed parameters
        precision = particleInputCI.precision
        user_particles_class = particleInputCI.user_particles_class
        particle_species = particleInputCI.particle_species

        self.position_coordinates = particleInputCI.position_coordinates
#        moved this out of ParticleDistributions_C:
#        self.position_coordinates = user_particles_class.position_coordinates
        self.dimension = len(self.position_coordinates)

        # These are the position coordinates at the start of a push
        initial_position_coordinates = [coord+'0' for coord in self.position_coordinates]

#        self.position_coordinates = particleInputCI.position_coordinates
        velocity_coordinates = ['u'+coord for coord in self.position_coordinates]
#        phase_coordinates = particleInputCI.position_coordinates + velocity_coordinates

        phase_coordinates = self.position_coordinates + initial_position_coordinates + velocity_coordinates

#        print 'particle_c... phase_coords = ', phase_coordinates
        
        # Count the species
        self.number_of_species = len(particle_species)
        if printFlag: print 'Particle_C: there are:', self.number_of_species, ' species'

        # Put the species names in a list called "names"
        self.species_names = [sp[0] for sp in particle_species]
        if printFlag: print 'Particle_C: species names:', self.species_names

        # Make reverse lookup dictionary index[] to give the species
        # index (starting from 1) given it's name.
        isp = 0
        self.index = {}
        for sn in self.species_names:
            isp += 1
            self.index[sn] = isp
            if printFlag: print 'Species', sn, 'is number', self.index[sn]

        # Put the user-defined plasma attributes in the following
        # dictionaries, which are indexed by the name of the species
        self.initial_distribution_type = {}
        self.initial_distribution_function = {}
        self.charge = {}
        self.dynamics = {}
        self.qom = {}

# don't know about this; not a fundamental number; just a particular
# initializer constant for a particular initialization

#        self.number_per_cell = {}
        self.explicit_species = []
        self.implicit_species = []
        self.neutral_species = []

        self.pseg_arr = {}
        self.particle_count = {}

        for sp in particle_species:
            species_name = sp[0]
            sp_dict = sp[1]

            # Process user input for initial particle distribution functions
            # key: 'initial_distribution_type'
            init_dist_type = sp_dict['initial_distribution_type']
            if init_dist_type:
                # Check for an 'attribute' in the
                # UserParticle/ParticleDistributions class that has
                # the same name as the species; that's the
                # user-specified initial distribution
#                if hasattr(UPD_M.UserParticleDistributions, init_dist):
                if hasattr(user_particles_class, species_name):
                    # store the name of the distribution function
                    self.initial_distribution_function[species_name] = getattr(user_particles_class, species_name)
                    if printFlag: print 'Particle_C: Initial distribution for', species_name, ' is the function of that name in ', user_particles_class
                # Write error message and exit if no distribution function exists
                else:
                    error_msg = "Particle_C: Need to define a particle distribution function %s in UserParticle.py for species %s " % (species_name, species_name)                                                                                                                                                                                                                                                                                                                       
                    sys.exit(error_msg)
            else:
                # There are no initial particles for this species
                    self.initial_distribution_function[species_name] = None

            self.initial_distribution_type[species_name] = init_dist_type

#            self.initial_distribution[species_name] = sp_dict['initial_distribution']

            # Process user input for fixed values of this species
            # key: 'charge'
            if printFlag: print 'Particle_C: sp_dict =', sp_dict
            self.charge[species_name] = sp_dict['charge']
            if printFlag: print 'Particle_C: charge for', species_name, 'is', self.charge[species_name]
            if sp_dict['charge'] == 0.0:
                self.neutral_species.append(species_name)

            # key: 'mass'
            self.qom[species_name] = sp_dict['charge']/sp_dict['mass']
            if printFlag: print 'Particle_C: charge-to-mass ratio for', species_name, 'is', self.qom[species_name]
            # key: 'number_per_cell'

# Should number_per_cell be here?  Eg if you have just some test particles
# Put it as part of the user's description of the distribution function
#            self.number_per_cell[species_name] = sp_dict['number_per_cell']
#            if echoFlag: print 'Particle_C: number per cell for ', species_name, ' is ', self.number_per_cell[species_name]

            # key: 'dynamics'
            self.dynamics[species_name] = sp_dict['dynamics']

            if sp_dict['dynamics'] == 'explicit':
                self.explicit_species.append(species_name)
            elif sp_dict['dynamics'] == 'implicit':
                self.implicit_species.append(species_name)
            else:
                error_msg = "Unknown type of dynamics ", sp_dict['dynamics'], ' for species ', sp_name

            # Process user input giving the particle-variable names and types
            # for each plasma species.  Allocate initial storage
            # for particles using segmented vectors indexed by the
            # species name.
            # key: 'attributes'
#            attribs = sp_dict['attributes']
            pvars = [coord for coord in phase_coordinates]
            pvars.append('weight')

#            pvartypes = [precision for coord in phase_coordinates]
            pvartypes = [precision for var in pvars]

            pvars.append('bitflags')
            pvartypes.append(np_M.int32)

            pvars.append('cell_index')
            pvartypes.append(np_M.int32) # The size determines how many local cells you can have.

            self.particle_dtype = {'names' : pvars, 'formats': pvartypes}
# just use ['bitflag']
#            self.bitflag_index = self.particle_dtype['names'].index('bitflags')

            if printFlag: print "Particle_C: particle metadata = %s" % self.particle_dtype

#            self.pseg_arr[species_name] = SA_M.SegmentedArray_C(segment_length, metadata)
            self.pseg_arr[species_name] = SA_M.SegmentedArray_C(self.SEGMENT_LENGTH, self.particle_dtype)

            # Initialize particle count for each species
            self.particle_count[species_name] = 0

        # May have particle trajectories
        self.trajCI = None

        # An scratch ndarray for one particle is used for trajectories
        self.one_particle_arr = np_M.empty(1, dtype=self.particle_dtype)

        # A scratch array that can hold: x,y,z, x0,y0,z0 (or subset)
        self.pcoords = np_M.empty(2*self.dimension, dtype=precision)
        # A scratch array that can hold: dx,dy,dz (or subset)
        self.dx = np_M.empty(self.dimension, dtype=precision)
        self.dx_in_cell = np_M.empty(self.dimension, dtype=precision)

        # Make a reusable array "self.negE" for computing -E at particle positions
        if particleInputCI.force_components is not None:
            Ecomps = particleInputCI.force_components
            force_precision = particleInputCI.force_precision
            Etypes = [force_precision for comp in Ecomps]
            Eseg_dict = {'names': Ecomps, 'formats': Etypes}
            self.negE = np_M.empty(self.SEGMENT_LENGTH, dtype=Eseg_dict)

        # Not used yet
        self.particle_integration_loop = particleInputCI.particle_integration_loop

        return

#class Particle_C(object):
    def initialize_distributions(self, printFlags):
        """Generate the initial particle distributions.
        """
#    pCI = runCI.particles # abbrev for particle Class Instance

        for sp in self.species_names:
            init_dist_type = self.initial_distribution_type[sp]
            if init_dist_type == 'listed':
                self.create_from_list(sp, printFlags[sp])
            elif self.initial_distribution_type == 'functional':
                self.create_from_functions(sp, printFlags[sp])
            elif self.initial_distribution_type == 'particle_file':
                self.create_from_file(sp, printFlags[sp])
            else:
                error_msg = "Unknown initial_distribution_type ", self.initial_distribution_type, " in Main for species ", sp
                sys.exit(error_msg)
        return
#    def initialize_distributions(self, printFlags):ENDDEF

#class Particle_C(object):
    def create_from_list(self, species_name, printFlag=False):
        """Generates particles for a species from a list provided by the user.
        """
        fncname = sys._getframe().f_code.co_name + '():'

        pseg_arrCI = self.pseg_arr[species_name] # The SegmentedArray_C object for this species
        
        # Set the current segmented array to be the 'out' array, since
        # we're going to write values to it.
#        pseg_arrCI.select_segmented_array('out')

        # The function that has the list of particles
        p_list_method = self.initial_distribution_function[species_name]

        number_of_macroparticles, particle_list = p_list_method(type='listed')

        # Check the length of the particle data by looking at the
        # datatype of the first segment of the segmented array.

        if len(particle_list[0]) != len(pseg_arrCI[0].dtype):
            error_msg = "Particle_Module.py: %s Species %s. Expecting particle data tuple %s, but data for first particle is: %s, which does not match up. Check UserParticles" %  (fncname, species_name, pseg_arrCI[0].dtype.names, particle_list[0])
#            print fncname, "Expect particle data for", pseg_arrCI[0].dtype.names
#            print "First particle is:", particle_list[0]
            sys.exit(error_msg)

        # Add the particles to storage and set up trajectories

#        print 'dtype = ', self.particle_dtype
#        print self.particle_dtype['names']
#        print 'index = ', self.particle_dtype['names'].index('bitflags')

        # put dynamics type here
        if self.trajCI is not None:
            if species_name in self.trajCI.explicit_species:
                dynamics_type = 'explicit'
            elif species_name in self.trajCI.implicit_species:
                dynamics_type = 'implicit'
            elif species_name in self.trajCI.neutral_species:
                dynamics_type = 'neutral'
            else:
                error_msg = "Particle_C:create_from_list: dynamics_type is unknown for species %s" % species_name
                sys.exit(error_msg)

        for i in range(number_of_macroparticles):
#            print 'species_name, particle_list[i] = ', species_name, particle_list[i]
            p, pindex = pseg_arrCI.put(particle_list[i])
            # Check if this particle has the trajectory flag turned on
#            if p[self.bitflag_index] & self.TRAJECTORY_FLAG == True:
            if p['bitflags'] & self.TRAJECTORY_FLAG != 0:
# or: p['bitflags'] should work?
                if self.trajCI is not None:
#                    print 'pindex for trajectory = ', pindex
                    self.trajCI.ParticleIdList[species_name].append(pindex)
                    self.trajCI.create_trajectory(species_name, dynamics_type)

#        if (printFlag): print fncname, "weight for ", species_name, " is ", weight
#        if (printFlag): print fncname, "bitflags for ", species_name, " is ", bitflags

        return
#    def create_from_list(self, species_name, printFlag=False):ENDDEF

#class Particle_C(object):
    def create_from_functions(self, species_name, printFlag = False):
        """Generates particles for a species from a list provided by the user.
        """
        fncname = sys._getframe().f_code.co_name + '():'

        pseg_arrCI = self.pseg_arrCI[species_name] # storage array for this species

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

        if (printFlag): print fncname, "weight for", species_name, "is", weight

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

#class Particle_C(object):
    def get_species_particle_count(self, species_name, printFlag = False):
        """Counts the particles for a given species.
        """
        fncname = sys._getframe().f_code.co_name + '():'

        pseg_arrCI = self.pseg_arr[species_name] # storage array for this species
        number_of_particles = pseg_arrCI.get_number_of_items()
        if printFlag: print fncname, species_name, 'has', number_of_particles, 'macroparticles'
# this should be done by the calling function?
#        self.particle_count[species_name] = number_of_particles

        return number_of_particles
#    def get_species_particle_count(self, species_name, printFlag = False):ENDDEF

#class Particle_C(object):
    def get_total_particle_count(self, printFlag = False):
        """Counts all the particles.
        """
        fncname = sys._getframe().f_code.co_name + '():'

        psum = 0
        for sn in self.species_names:
            pseg_arrCI = self.pseg_arr[sn] # storage array for this species
            npar = pseg_arrCI.get_number_of_items()
            psum += npar

        if printFlag: print fncname, 'Total number of macroparticles in', len(self.species_names), 'species is', psum
        return psum

#    def get_total_particle_count(self, printFlag = False):ENDDEF

#class Particle_C(object):
    def compute_mesh_cell_indices(self, meshCI):
        """Compute the cell index for each particle.

           :param meshCI: A Mesh_C object containing the mesh attributes.

           :returns: None

        """

        for sn in self.species_names:
            psaCI = self.pseg_arr[sn] # segmented array for this species
            (np_seg, pseg) = psaCI.init_out_loop()
#            pseg = psaCI.get_next_segment()

            while isinstance(pseg, np_M.ndarray):
#                for ip in xrange(pseg.size):
                for ip in xrange(np_seg):
                    pseg[ip]['cell_index'] = meshCI.compute_cell_index(pseg[ip])
#                    print 'ip, index =', ip, pseg[ip]['cell_index']
                (np_seg, pseg) = psaCI.get_next_segment('out')

        return
#    def compute_mesh_cell_indices(self, meshCI):ENDDEF

#class Particle_C(object):
    def move_particles_in_electrostatic_field(self, dt, neg_electric_field):
        """Apply electric force to particles.  Compute the change in
           velocity and position in time dt. Use an explicit
           method to integrate the orbit.

           Arguments:
               dt: time interval.

               neg_electric_field: a Field_C object

        """

        fncname = sys._getframe().f_code.co_name + '():'

        meshCI = neg_electric_field.meshCI

        # Scratch space
        pcoords = self.pcoords # x,y,z, x0,y0,z0 (or subset)
        dx = self.dx # dx is the distance moved in one step.
#        p_arr = self.one_particle_arr[0]

        dim = self.dimension

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
#            (np_seg, pseg_in, pseg_out) = psaCI.init_inout_loop()
            (np_seg, pseg_in) = psaCI.init_inout_loop()
#            psaCI.init_segment_loop()
#            pseg = psaCI.get_next_segment()
#            while pseg != None:

            particle_count = 0
            # ip_out counts through the 'out' array
            ip_out = 0
            while isinstance(pseg_in, np_M.ndarray):
#            print pseg['z']
#                print 'particles_mod: particles module: pseg = ', pseg
                # Compute electric field for each particle

#                print "Before:"
#                print "position:", pseg['x'], pseg['y'] #, pseg['z']
#                print "velocity:", pseg['ux'], pseg['uy'] #, pseg['uz']

                # This is the negative electric field
                neg_electric_field.interpolate_field_to_points(pseg_in, self.negE)
                
                # Truncate negE to the number of particles to get the
                # += operations below to work: These operations need
                # the vectors to be the same length.
#                Eseg = self.negE[0:pseg_in.size]
                Eseg = self.negE[0:np_seg]
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

                for ip_in in xrange(np_seg):
#                for ip in xrange(pseg.size):

                    # pseg[i] is 'x', 'y', 'z', 'x0',..., 'ux', 'uy',... values of ith item
                    # So pseg[i][0:3] is 'x', 'y', 'z'.
                    # Can't use slice syntax here, because the array data are not homogeneous.

                    # skip deleted particles
                    if pseg_in[ip_in]['bitflags'] & self.DELETE_FLAG != 0: continue
                    
                    # Only gets the next 'out' segment if there are more 'in' particles to do
                    if ip_out == 0: pseg_out = psaCI.get_next_out_segment()

                    # Copy the 'in' particle to a tmp
#                    p_arr = pseg_in[ip_in]
#                    print 'p_arr = ', p_arr
                    pseg_out[ip_out] = pseg_in[ip_in] # This copies everything over, to ensure weights, flags, etc. get copied.
#                    print 'pseg_out = ', pseg_out[ip_out]
#                    print 'same?', pseg_out[ip_out] is pseg_in[ip_in]


                    # Accelerate the particle with the field components in Eseg
                    for comp in Eseg.dtype.names:
                        ucomp = 'u'+comp
#                        pseg_out[ip_out][ucomp] = pseg_in[ip_in][ucomp] + qmdt*Eseg[ip_in][comp]
                        pseg_out[ip_out][ucomp] += qmdt*Eseg[ip_in][comp] # Index OK?

# Instead of this, could do ifs on the dimension of Eseg
                    """                    
                    pseg_out[ip_out]['ux'] = pseg_in[ip_in]['ux'] + qmdt*Eseg[ip_in]['x']
                    pseg_out[ip_out]['uy'] = pseg_in[ip_in]['uy'] + qmdt*Eseg[ip_in]['y']
                    pseg_out[ip_out]['uz'] = pseg_in[ip_in]['uz'] + qmdt*Eseg[ip_in]['z']
                    """                    

# Need the dim value

                    # Move the particle
                    # NON-RELATIVISTIC: v and u are the same

                    for coord in self.position_coordinates:
                        coord0 = coord+'0'
#                        pseg_out[ip_out][coord0] = pseg_in[ip_in][coord] # Save the starting positions
                        pseg_out[ip_out][coord0] = pseg_out[ip_out][coord] # Save the starting positions
                        ucomp = 'u'+coord
#                        pseg_out[ip_out][coord] = pseg_in[ip_in][coord] + pseg_in[ip_in][ucomp]*dt
                        pseg_out[ip_out][coord] += pseg_out[ip_out][ucomp]*dt

                    """                    
                    pseg_out[ip_out]['x0'] = pseg_in[ip_in]['x']
                    pseg_out[ip_out]['y0'] = pseg_in[ip_in]['y']
                    pseg_out[ip_out]['z0'] = pseg_in[ip_in]['z']

                    pseg_out[ip_out]['x'] = pseg_in[ip_in]['x'] + pseg_out[ip_out]['ux']*dt
                    pseg_out[ip_out]['y'] = pseg_in[ip_in]['y'] + pseg_out[ip_out]['uy']*dt
                    pseg_out[ip_out]['z'] = pseg_in[ip_in]['z'] + pseg_out[ip_out]['uz']*dt
                    """

                    # Check if the particle is still on the meshed region

#                    print 'ip, index =', ip, pseg[ip]['cell_index']
                    p_cell_index = pseg_out[ip_out]['cell_index']
#                    print fncname, ": ip, pindex", ip, p_cell_index, "cell index:", meshCI.compute_cell_index(pseg[ip])
                    if meshCI.is_inside(pseg_out[ip_out], p_cell_index):
                        pass
                    else:
                        # The particle has left this cell.  We may
                        # need to track it across each facet in case
                        # there's a boundary-condition on that facet.

                        # Loop on facets of this cell to find which
                        # one it crossed.

                        # The particle may have hit a boundary

                        # To detect which cell facet it went through:

                        # 1. Compute the distance a particle moves
                        # in the directions of a facet (Lf).

                        # 2. If it moved closer to it, compute the
                        # distance to the facet (Df).

                        # 3. If Df < Lf, the particle crossed the
                        # plane that the facet is in.

                        # 4. For the crossed planes, the facet
                        # crossed is the one with the smallest
                        # ratio Df/Lf, since it reached that facet
                        # first.

#                        pcoords[:] = pseg_out[ip_out] # Alias to the position coordinates
                        '''
                        i=0
                        for coord in self.position_coordinates:
                            pcoords[i] = pseg_out[ip_out][coord]
                            coord0 = coord+'0'
                            pcoords[i+dim] = pseg_out[ip_out][coord0]
                            i+=1
                        '''
                        for i in range(2*dim):
                            pcoords[i] = pseg_out[ip_out][i]

                        dx = pcoords[0:dim] - pcoords[dim:2*dim] # Move vector

                        found_facet = False

#                        facet = meshCI.find_facet(pcoords, pcoords[dim+1:2*dim)
                        facet = meshCI.find_facet(pcoords, dx, p_cell_index)

# See save10/Particle_Module.py for a different search

#                        if found_cell != True: print "Particle is not in nearby cells; using BB search"
#                        ci = meshCI.compute_cell_index(pseg[ip])
#                        print "Found particle is in cell", ci        
#                        pseg[ip]['cell_index'] = ci

                    # Don't need this since we just look up the cell
                    # index when computing negE above
                    # pseg[ip]['cell_index'] = compute_cell_index(Point(p))

                    ip_out += 1
                    # Check if we've reached the end of this segment.  If
                    # so, we need to start writing on a new segment.
                    if (ip_out == self.SEGMENT_LENGTH):
                        particle_count += ip_out
                        ip_out = 0 # This will cause get_next_out_segment() to be called
                                   # above if there are more 'in' particles to be processed.

                # Done with this segment.
                # Get the next one, if it exists.
                (np_seg, pseg_in) = psaCI.get_next_segment('in')
#                pseg = psaCI.get_next_segment()
                # Loop over particles in this segment ends

            # Set values that will be used to keep track of the particle arrays
            particle_count += ip_out
            psaCI.set_number_of_items('out', particle_count)
            # Loop over segmented array ends

# Compute new density here?

        # Loop over Species ends
        return
#    def move_particles_in_electrostatic_field(self, dt, neg_electric_field):ENDDEF

#class Particle_C(object):
    def move_neutral_particles(self, dt, meshCI):
        """Move neutral particles on a mesh.
           
           Compute change in position in time dt. Use an explicit
           method to integrate the orbit.

           Arguments:
               dt: time interval.

               meshCI: a Mesh_C object

        """

        fncname = sys._getframe().f_code.co_name + '():'

        # Scratch space
        pcoords = self.pcoords # x,y,z, x0,y0,z0 (or subset)
        dx = self.dx # dx is the distance moved in one step.
#        p_arr = self.one_particle_arr[0]

        dim = self.dimension

        for sn in self.neutral_species:

            # Invariant parameters
#            print 'sn = ', sn
        
#        nlocal = self.particles.get_species_particle_count(sp, echoFlag = False)
            psaCI = self.pseg_arr[sn] # segmented array for this species

            if self.get_species_particle_count(sn) == 0: continue

#           Move all the particles in this species
            (np_seg, pseg_in) = psaCI.init_inout_loop()
            particle_count = 0
            # ip_out counts through the 'out' array
            ip_out = 0
            while isinstance(pseg_in, np_M.ndarray):
                for ip_in in xrange(np_seg):
                    # pseg[i] is 'x', 'y', 'z', 'x0',..., 'ux', 'uy',... values of ith item
                    # So pseg[i][0:3] is 'x', 'y', 'z'.
                    # Can't use slice syntax here, because the array data are not homogeneous.

                    # skip deleted particles
                    if pseg_in[ip_in]['bitflags'] & self.DELETE_FLAG != 0: continue
                    
                    # Only gets the next 'out' segment if there are more 'in' particles to do
                    if ip_out == 0: pseg_out = psaCI.get_next_out_segment()

                    pseg_out[ip_out] = pseg_in[ip_in] # This copies everything over, to ensure weights, flags, etc. get copied.

                    # Move the particle
                    # NON-RELATIVISTIC: v and u are the same

                    for coord in self.position_coordinates:
                        coord0 = coord+'0'
                        pseg_out[ip_out][coord0] = pseg_out[ip_out][coord] # Save the starting positions
                        ucomp = 'u'+coord
                        pseg_out[ip_out][coord] += pseg_out[ip_out][ucomp]*dt
                        
#                    print "x0=", pseg_out[ip_out]['x0'],"x=", pseg_out[ip_out]['x']
#                    print "pseg_out", pseg_out[ip_out]

                    # Check if the particle is still on the meshed region

#                    print 'ip, index =', ip, pseg[ip]['cell_index']
                    p_cell_index = pseg_out[ip_out]['cell_index']
#                    print fncname, ": ip, pindex", ip, p_cell_index, "cell index:", meshCI.compute_cell_index(pseg[ip])
                    if meshCI.is_inside(pseg_out[ip_out], p_cell_index):
                        pass
                    else:
                        # The particle has left this cell.  We may
                        # need to track it across each facet in case
                        # there's a boundary-condition on that facet.
#                        print fncname, "particle has migrated"
                        i=0
                        for coord in self.position_coordinates:
                            pcoords[i] = pseg_out[ip_out][coord]
                            coord0 = coord+'0'
                            pcoords[i+dim] = pseg_out[ip_out][coord0]
                            i+=1
#                        pcoords[:] = pseg_out[ip_out] # Alias to the position coordinates
                        dx = pcoords[0:dim] - pcoords[dim:2*dim] # Move vector

#                        facet = meshCI.find_facet(pcoords, pcoords[dim+1:2*dim)
#                        print fncname, "pcoords=", pcoords, "dx=", dx, "p_cell_index=", p_cell_index

                        (facet, path_fraction) = meshCI.find_facet(pcoords, dx, p_cell_index)
                        if facet != meshCI.NO_FACET:
                            # Update the cell index to the new cell
# NB: This may not be the final cell. FIXME, with a while loop???
                            pseg_out[ip_out]['cell_index'] = meshCI.cell_neighbor_dict[p_cell_index][facet]
#                            print "cell_index updated to", pseg_out[ip_out]['cell_index']
                        else:
                            print "facet crossed is", facet

                    # Don't need this since we just look up the cell
                    # index when computing negE above
                    # pseg[ip]['cell_index'] = compute_cell_index(Point(p))

                    ip_out += 1
                    # Check if we've reached the end of this segment.  If
                    # so, we need to start writing on a new segment.
                    if (ip_out == self.SEGMENT_LENGTH):
                        particle_count += ip_out
                        ip_out = 0 # This will cause get_next_out_segment() to be called
                                   # above if there are more 'in' particles to be processed.

                # Done with this segment.
                # Get the next one, if it exists.
                (np_seg, pseg_in) = psaCI.get_next_segment('in')
#                pseg = psaCI.get_next_segment()
                # Loop over particles in this segment ends

            # Set values that will be used to keep track of the particle arrays
            particle_count += ip_out
            psaCI.set_number_of_items('out', particle_count)
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
            while isinstance(pseg, np_M.ndarray):
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
    def move_particles_in_uniform_fields(self, species_name, ctrlCI, printFlag = False):
        """Apply electric field ctrlCI.E0 to particles.  Compute the
           resulting change in particle velocities and positions in time
           ctrlCI.dt.
        """
# Use numpy array syntax to do this, instead of loops; see Susie; See HPL Sec 4.2.2

        dt = ctrlCI.dt
        E0 = ctrlCI.E0

        # Invariant parameters
        
        qmdt = self.qom[species_name]*dt

#        nlocal = self.particles.get_species_particle_count(sp, echoFlag = False)

        # Accelerate the particles
        psaCI = self.pseg_arr[species_name] # segmented array for this species

        # This loop moves particles, so swap the in/out particle arrays
#        (np_seg, pseg_in, pseg_out) = psaCI.init_inout_loop()
        (np_seg, pseg_in) = psaCI.init_inout_loop()

#        print "move_particles_in_uniform_fields: np_seg, pseg_in, pseg_out:", np_seg, pseg_in, pseg_out
        # Get the first segment of particles to move
#        pseg = psaCI.get_next_segment()

        # Loop on segments until a None value is returned
        particle_count = 0
        # ip_out counts through the 'out' array
        ip_out = 0
        while isinstance(pseg_in, np_M.ndarray):
#            print "type pseg = ", type(pseg)
#            print pseg['z']
#            print pseg
            # Accelerate this block of particles

            # Loop on the particles in this segment.
            # If a particle has been deleted, skip it.

            # ip_in counts through the 'in' segment
            for ip_in in xrange(np_seg):

                # Skip deleted particles
                if pseg_in[ip_in]['bitflags'] & self.DELETE_FLAG != 0: continue

                # Only gets the next 'out' segment if there are more 'in' particles to do
                if ip_out == 0: pseg_out = psaCI.get_next_out_segment()

                pseg_out[ip_out] = pseg_in[ip_in] # This copies everything over, to ensure weights, flags, etc. get copied.

                # Accelerate the particle
                pseg_out[ip_out]['ux'] = pseg_in[ip_in]['ux'] + qmdt*E0.x
                pseg_out[ip_out]['uy'] = pseg_in[ip_in]['uy'] + qmdt*E0.y
                pseg_out[ip_out]['uz'] = pseg_in[ip_in]['uz'] + qmdt*E0.z

                # Move the particle
                # NON-RELATIVISTIC: v and u are the same
                pseg_out[ip_out]['x'] = pseg_in[ip_in]['x'] + pseg_out[ip_out]['ux']*dt
                pseg_out[ip_out]['y'] = pseg_in[ip_in]['y'] + pseg_out[ip_out]['uy']*dt
                pseg_out[ip_out]['z'] = pseg_in[ip_in]['z'] + pseg_out[ip_out]['uz']*dt

                ip_out += 1
                # Check if we've reached the end of this segment.  If
                # so, we need to start writing on a new segment.  But
                # don't increment to a new 'out' segment unless there
                # are actually more 'in' particles to process.
                if (ip_out == self.SEGMENT_LENGTH):
                    particle_count += ip_out
                    ip_out = 0 # This will cause get_next_out_segment() to be called
                               # above if there are more 'in' particles to be processed.

    #            print pseg['z']
            (np_seg, pseg_in) = psaCI.get_next_segment('in')

        # Set values that will be used to keep track of the particle arrays
#        print "Species", species_name, "np_seg=", np_seg, "pseg_in", pseg_in, "particle_count", particle_count
        particle_count += ip_out
        psaCI.set_number_of_items('out', particle_count)

        return
#    def move_particles_in_uniform_fields(self, species_name, ctrlCI, printFlag = False):ENDDEF

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
#    def move_particles_in_uniform_fields(self, species_name, ctrlCI, printFlag = False):ENDDEF


#class Particle_C(object):
    def record_trajectory_data(self, neg_electric_field):
        """Save the requested particle trajectory data.
        """

        trajCI = self.trajCI
        p_arr = self.one_particle_arr
#        print 'p_arr.dtype.names = ', p_arr.dtype.names
        count = self.trajCI.count

        # Loop on trajectories of explicit particles
        for sp in trajCI.explicit_species:
            pseg_arrCI = self.pseg_arr[sp] # The SA for this species
            for i in xrange(len(trajCI.ParticleIdList[sp])):
#            for ip in trajCI.ParticleIdList[sp]:
                ip = trajCI.ParticleIdList[sp][i]
                # Retrieve particle using its local index
                p_arr[0] = pseg_arrCI.get(ip)
#                print "record_trajectory_data: p_arr[0] = ", p_arr[0]

                # Compute the force on this particle

#                E_arr = field_particleCI.compute_E_at_particles(p_arr)

                neg_electric_field.interpolate_field_to_points(p_arr, self.negE)

                E_arr = self.negE[0:p_arr.size]
                # Flip the sign to get correct E. (This is a vector
                # operation on each component of E.)
                for n in E_arr.dtype.names:
                    E_arr[n] *= -1.0

#                print "E_arr.dtype.names =", E_arr.dtype.names

                # Copy the particle values into the trajectory
                for comp in trajCI.explicitDict['names']:
#                    print 'comp = ', comp
                    if comp in p_arr.dtype.names:
#                        trajCI.DataList[sp][i][count][comp] = p_arr[0][comp]
                        trajCI.DataList[sp][i][comp][count] = p_arr[0][comp]
                    if comp in E_arr.dtype.names:
                        Ecomp = 'E'+comp
                        trajCI.DataList[sp][i][Ecomp][count] = E_arr[0][comp]
#                        trajCI.DataList[sp][i][count][Ecomp] = E_arr[0][comp]

#        print 'trajCI.DataList: ', trajCI.DataList['testelectrons'][0]

        # Loop on trajectories of implicit particles
        
        self.trajCI.count += 1

        return
#    def record_trajectory_data(self, neg_electric_field):ENDDEF

#class Particle_C(object):
    def compute_charge_density(self, fCI):
        """Apply electric force to particles.  Compute the change in
           velocity and position in time dt. Use an explicit method to
           integrate the orbit.

           Arguments:
               dt: time interval.

               fpCI: an object that can provide method to compute
               density from discrete particles.
               
        """

        isp = 0
        for sn in self.species_names:

            # Invariant parameters
#            print 'sn = ', sn
            charge = self.charge[sn]
        

# charge on the particle:
# compute density first, then accumulate into a charge-density

# species number?


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
            while isinstance(pseg, np_M.ndarray):
                fCI.compute_rho_from_particles(isp, charge, pseg)

            isp += 1

            return
#    def compute_charge_density(self, fCI):ENDDEF

#class Particle_C(object):ENDCLASS
