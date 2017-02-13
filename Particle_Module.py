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

#    def __init__(self, species_input, phase_coordinates, precision, segment_length, user_particles_class, echoFlag):


#    def __init__(self, species_input, phase_coordinates, precision, user_particles_class, echoFlag):
#class Particle_C(object):
    def __init__(self, particleInputCI, printFlag = False):
        """Take a list of kinetic species provided by the user and create the initial plasma
        """
        # Set local variables from passed parameters
        precision = particleInputCI.precision
        user_particles_class = particleInputCI.user_particles_class
        particleSpecies = particleInputCI.particle_species

        # The spatial coordinates of a particle
        self.position_coordinates = particleInputCI.position_coordinates
        self.particle_dimension = len(self.position_coordinates)

        # These are the position coordinates at the start of a push
        initial_position_coordinates = [coord+'0' for coord in self.position_coordinates]

#        self.position_coordinates = particleInputCI.position_coordinates
        velocity_coordinates = ['u'+coord for coord in self.position_coordinates]
#        phase_coordinates = particleInputCI.position_coordinates + velocity_coordinates

        phase_coordinates = self.position_coordinates + initial_position_coordinates + velocity_coordinates
#        print 'particle_c... phase_coords = ', phase_coordinates

        # This is for a reference to a UserMesh_C object for particles
        self.pmeshCI = None
        
        # This is for a reference to a ParticleMeshBoundaryConditions_C
        # object that handles particle boundary conditions.
        self.pmesh_bcCI = None

        # Count the species
        self.number_of_species = len(particleSpecies)
        if printFlag: print 'Particle_C: there are:', self.number_of_species, ' species'

        # Put the species names in a list called "names"
        self.species_names = [sp[0] for sp in particleSpecies]
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

        for sp in particleSpecies:
            species_name = sp[0]
            sp_dict = sp[1]

            # Process user input for initial particle distribution functions
            # key: 'initial_distribution_type'
            init_dist_type = sp_dict['initial_distribution_type']
            if init_dist_type is not None:
                # Look in user_particles_class for a function
                # ('attribute') that has the same name as the particle
                # species; that's the user-specified initial
                # distribution.
#                if hasattr(UPD_M.UserParticleDistributions, init_dist):
                if hasattr(user_particles_class, species_name):
                    # store the name of the distribution function
                    self.initial_distribution_function[species_name] = getattr(user_particles_class, species_name)
                    if printFlag: print 'Particle_C: Initial distribution for', species_name, ' is the function of that name in ', user_particles_class
                # Write error message and exit if no distribution function exists
                else:
                    errorMsg = "Particle_C: Need to define a particle distribution function %s in UserParticle.py for species %s " % (species_name, species_name)
                    sys.exit(errorMsg)
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
                errorMsg = "Unknown type of dynamics ", sp_dict['dynamics'], ' for species ', sp_name

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

        # This is for a reference to a Trajectory_C object to handle
        # particles that have the TRAJECTORY_FLAG bit turned on.  It
        # needs to have that object before such particles are
        # encountered.
        self.trajCI = None

        # An scratch ndarray for one particle is used for trajectories
        self.one_particle_arr = np_M.empty(1, dtype=self.particle_dtype)

        # A scratch array that can hold: x,y,z, x0,y0,z0 (or subset)
        self.pcoords = np_M.empty(2*self.particle_dimension, dtype=precision)
        # A scratch array that can hold: dx,dy,dz (or subset)
        self.dx = np_M.empty(self.particle_dimension, dtype=precision)
        self.dx_in_cell = np_M.empty(self.particle_dimension, dtype=precision)

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

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

        for sp in self.species_names:
            init_dist_type = self.initial_distribution_type[sp]
            if init_dist_type == 'listed':
                self.create_from_list(sp, printFlags[sp])
            elif self.initial_distribution_type == 'functional':
                self.create_from_functions(sp, printFlags[sp])
            elif self.initial_distribution_type == 'particle_file':
                self.create_from_file(sp, printFlags[sp])
            else:
                errorMsg = "Unknown initial_distribution_type ", self.initial_distribution_type, " in Main for species ", sp
                sys.exit(errorMsg)
        return
#    def initialize_distributions(self, printFlags):ENDDEF

#class Particle_C(object):
    def create_from_list(self, species_name, printFlag=False):
        """Generates particles for a species from a list provided by the user.
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

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
            errorMsg = "Particle_Module.py: %s Species %s. Expecting particle data tuple %s, but data for first particle is: %s, which does not match up. Check UserParticles" %  (fncName, species_name, pseg_arrCI[0].dtype.names, particle_list[0])
#            print fncName, "Expect particle data for", pseg_arrCI[0].dtype.names
#            print "First particle is:", particle_list[0]
            sys.exit(errorMsg)

        # Add the particles to storage and set up trajectories

#        print 'dtype = ', self.particle_dtype
#        print self.particle_dtype['names']
#        print 'index = ', self.particle_dtype['names'].index('bitflags')

        # The dynamics type determines what trajectory variables are available
        if self.trajCI is not None:
            if species_name in self.trajCI.explicit_species:
                dynamics_type = 'explicit'
            elif species_name in self.trajCI.implicit_species:
                dynamics_type = 'implicit'
            elif species_name in self.trajCI.neutral_species:
                dynamics_type = 'neutral'
            else:
                errorMsg = "Particle_C:create_from_list: dynamics_type is unknown for species %s" % species_name
                sys.exit(errorMsg)

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
                else:
# Instead of printing this message, a trajCI object could be created here.
                    print fncName, "*** DT Warning: A trajectory flag is on, but no trajectory object has been created yet. ***"

#        if (printFlag): print fncName, "weight for ", species_name, " is ", weight
#        if (printFlag): print fncName, "bitflags for ", species_name, " is ", bitflags

        return
#    def create_from_list(self, species_name, printFlag=False):ENDDEF

#class Particle_C(object):
    def create_from_functions(self, species_name, printFlag = False):
        """Generates particles for a species from a list provided by the user.
        """
        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

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

        if (printFlag): print fncName, "weight for", species_name, "is", weight

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
        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

        pseg_arrCI = self.pseg_arr[species_name] # storage array for this species
        number_of_particles = pseg_arrCI.get_number_of_items()
        if printFlag: print fncName, species_name, 'has', number_of_particles, 'macroparticles'
# this should be done by the calling function?
#        self.particle_count[species_name] = number_of_particles

        return number_of_particles
#    def get_species_particle_count(self, species_name, printFlag = False):ENDDEF

#class Particle_C(object):
    def get_total_particle_count(self, printFlag = False):
        """Counts all the particles.
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

        psum = 0
        for sn in self.species_names:
            pseg_arrCI = self.pseg_arr[sn] # storage array for this species
            npar = pseg_arrCI.get_number_of_items()
            psum += npar

        if printFlag: print fncName, 'Total number of macroparticles in', len(self.species_names), 'species is', psum
        return psum

#    def get_total_particle_count(self, printFlag = False):ENDDEF

#class Particle_C(object):
    def compute_mesh_cell_indices(self):
        """Compute the cell index for each particle.

           :returns: None

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

        if self.pmeshCI is not None:
            for sn in self.species_names:
                psaCI = self.pseg_arr[sn] # segmented array for this species
                (npSeg, pseg) = psaCI.init_out_loop()

                while isinstance(pseg, np_M.ndarray):
    #                for ip in xrange(pseg.size):
                    for ip in xrange(npSeg):
                        pseg[ip]['cell_index'] = self.pmeshCI.compute_cell_index(pseg[ip])
    #                    print 'ip, index =', ip, pseg[ip]['cell_index']
# Check that is_inside() confirms the cell index:
                        if not self.pmeshCI.is_inside(pseg[ip], pseg[ip]['cell_index']):
                            errorMsg = "%s !!! is_inside() check failed for particle %d" % (fncName, ip)
                            sys.exit(errorMsg)
#                        else:
#                            print fncName, "*** is_inside check passes for particle", pseg[ip], "***"
                    (npSeg, pseg) = psaCI.get_next_segment('out')
        else:
            print fncName, "!!! The reference to pmeshCI is None!!!"

        return
#    def compute_mesh_cell_indices(self):ENDDEF

#class Particle_C(object):
    def move_particles_in_electrostatic_field(self, dt, neg_electric_field):
        """Apply electric force to particles.  Compute the change in
           velocity and position in time dt. Use an explicit
           method to integrate the orbit.

           Arguments:
               dt: time interval.

               neg_electric_field: a Field_C object

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

#        meshCI = neg_electric_field.meshCI
        pmeshCI = self.pmeshCI

        # Scratch space
        pCoords = self.pcoords # x,y,z, x0,y0,z0 (or subset)
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
            (npSeg, psegIn) = psaCI.init_inout_loop()
            particleCount = 0
            # ipOut counts through the 'out' array
            ipOut = 0
            while isinstance(psegIn, np_M.ndarray):
#            print pseg['z']
#                print 'particles_mod: particles module: pseg = ', pseg
                # Compute electric field for each particle

#                print "Before:"
#                print "position:", pseg['x'], pseg['y'] #, pseg['z']
#                print "velocity:", pseg['ux'], pseg['uy'] #, pseg['uz']

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
                    if psegIn[ipIn]['bitflags'] & self.DELETE_FLAG != 0: continue
                    
                    # Only gets the next 'out' segment if there are more 'in' particles to do
                    if ipOut == 0: psegOut = psaCI.get_next_out_segment()

                    # Copy the 'in' particle to a tmp
#                    p_arr = psegIn[ipIn]
#                    print 'p_arr = ', p_arr
                    psegOut[ipOut] = psegIn[ipIn] # This copies everything over, to ensure weights, flags, etc. get copied.
#                    print 'psegOut = ', psegOut[ipOut]
#                    print 'same?', psegOut[ipOut] is psegIn[ipIn]


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

#                    print fncName, ": ip, pindex", ip, pCellIndex, "cell index:", pmeshCI.compute_cell_index(pseg[ip])
                    mLastFacet = pmeshCI.NO_FACET
                    facetCrossCount = 0
                    while not pmeshCI.is_inside(psegOut[ipOut], pCellIndex):
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
                            pCoords[i] = psegOut[ipOut][coord] # Present position
                            coord0 = coord+'0'
                            pCoords[i+pDim] = psegOut[ipOut][coord0] # Start of path in this cell
                            i+=1
                        '''
                        # Replace above by:
                        for i in range(2*pDim):
                            pCoords[i] = psegOut[ipOut][i]

                        dx = pCoords[0:pDim] - pCoords[pDim:2*pDim] # Move vector

# See save10/Particle_Module.py for a different search, where only nearby cells are looked at.
                        (cFacet, pathFraction) = pmeshCI.find_facet(pCoords[pDim:2*pDim], dx, pCellIndex)
                        if cFacet != pmeshCI.NO_FACET:

#                        if found_cell != True: print "Particle is not in nearby cells; using BB search"
#                        ci = pmeshCI.compute_cell_index(pseg[ip])
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
                            mFacet = pmeshCI.cell_entity_index_dict['facet'][pCellIndex][cFacet]
# mFacet should never be the same as the last facet crossed: check this
                            if mFacet == mLastFacet:
                                errorMsg = "%s The mesh index of the facet crossed is %d, the same as the last facet crossed. This should not happen since the particle cannot cross the same facet twice in one move!" % (fncName, mFacet)
                                sys.exit(errorMsg)
                            else:
                                mLastFacet = mFacet
                            # ...and get the value of the facet marker.
# mFacet is a numpy.uint32, but the FacetFunction wants an int argument.
#                            print "type is:", type(mFacet)
                            facValue = pmeshCI.particle_boundary_marker[mFacet]
                            # Check if this facet has a non-zero marker
                            if facValue != 0:
                                # Convert the int marker to its string value...
                                bFlag = str(facValue)
#                                print "bFlag is", bFlag
                                # ...and call the function associated with this value.

# Maybe this needs the whole particle object?  It could be the end-of-the line for the particle.
# Need to check what needs to be in the following arg list:
                                self.pmesh_bcCI.bcFunctionDict[bFlag][sn](psegOut[ipOut], sn, mFacet)

# Not sure if this else is needed:                            else:

                            # Update the cell index to the new cell
                            pCellIndex = pmeshCI.cell_neighbor_dict[pCellIndex][cFacet]
                            psegOut[ipOut]['cell_index'] = pCellIndex
                            # If the particle has left the grid, exit this 'while' loop.
                            if pCellIndex == pmeshCI.NO_CELL: break # Exit the 'while' loop
                            # If pathFraction is zero, this is a flag
                            # to verify that the index is correct
#                            if pathFraction == 0.0:
                        else:
                            errorMsg = "%s The cell index of the facet crossed is %d. This should not happen since the particle has left its initial cell cell!" % (fncName, cFacet)
                            sys.exit(errorMsg)
#                       END:if cFacet != pmeshCI.NO_FACET:
#
#                   END:while not pmeshCI.is_inside(psegOut[ipOut], pCellIndex)

                    # Don't need this since we just look up the cell
                    # index when computing negE above
                    # pseg[ip]['cell_index'] = compute_cell_index(Point(p))

                    # Check that this particle has not been deleted before
                    # incrementing the 'out' particle counter
                    if psegOut[ipOut]['bitflags'] & self.DELETE_FLAG == 0:
                        ipOut += 1
                    else:
                        print "Particle", ipIn, "of species", sn, "has been deleted."

                    # Check if we've reached the end of this segment.  If
                    # so, we need to start writing on a new segment.
                    if (ipOut == self.SEGMENT_LENGTH):
                        particleCount += self.SEGMENT_LENGTH
                        ipOut = 0 # This will cause get_next_out_segment() to be called
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
           
           Compute change in position in time dt. Use an explicit
           method to integrate the orbit.

           Arguments:
               dt: time interval.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

        # Scratch space
        pCoords = self.pcoords # x,y,z, x0,y0,z0 (or subset)
        dx = self.dx # dx is the distance moved in one step.
#        p_arr = self.one_particle_arr[0]

        pmeshCI = self.pmeshCI
        pDim = self.particle_dimension

        for sn in self.neutral_species:

            # Invariant parameters
#            print 'sn = ', sn
        
#        nlocal = self.particles.get_species_particle_count(sp, echoFlag = False)
            psaCI = self.pseg_arr[sn] # segmented array for this species

            if self.get_species_particle_count(sn) == 0: continue

#           Move all the particles in this species
            (npSeg, psegIn) = psaCI.init_inout_loop()
            particleCount = 0
            # ipOut counts particles being written to the current
            # 'out' segment.
            ipOut = 0
            while isinstance(psegIn, np_M.ndarray):
                for ipIn in xrange(npSeg):
                    # pseg[i] is 'x', 'y', 'z', 'x0',..., 'ux', 'uy',... values of ith item
                    # So pseg[i][0:3] is 'x', 'y', 'z'.
                    # Can't use slice syntax here, because the array data are not homogeneous.

                    # skip deleted particles
                    if psegIn[ipIn]['bitflags'] & self.DELETE_FLAG != 0: continue
                    
                    # Only gets the next 'out' segment if there are more 'in' particles to do
                    if ipOut == 0: psegOut = psaCI.get_next_out_segment()

# Check if the following copies, or resets the reference. Ans: it COPIES.
                    psegOut[ipOut] = psegIn[ipIn]
# Also, check if the following copies, or resets the reference. Ans: it COPIES.
#                    p_out = psegOut[ipOut]
# To get a REFERENCE, need to use SLICING syntax:
#                    p_out = psegIn[ipIn][:] # This copies everything over, to ensure weights, flags, etc. get copied.
# But you can't get a reference to a single element of psegIn[]. See:
# http://stackoverflow.com/questions/23654088/reference-of-a-single-numpy-array-element
# See also HPL p. 137. b=a[1,:], but b is still an array.
# What about getting a ref to one particle returned from a function? That's possible because it uses the stack?

                    # Move the particle
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
#                    print fncName, ": ip, pindex", ipOut, pCellIndex, "cell index:", pmeshCI.compute_cell_index(psegOut[ipOut])

                    # Loop until the particle is in the current cell
                    mLastFacet = pmeshCI.NO_FACET
                    facetCrossCount = 0
                    while not pmeshCI.is_inside(psegOut[ipOut], pCellIndex):
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
                            pCoords[i] = psegOut[ipOut][coord] # Present position
                            coord0 = coord+'0'
                            pCoords[i+pDim] = psegOut[ipOut][coord0] # Start of path in this cell
                            i+=1
#                        pCoords[:] = psegOut[ipOut] # Alias to the position coordinates
                        dx = pCoords[0:pDim] - pCoords[pDim:2*pDim] # Move vector

# could return the crossing-point in pCoords[]
                        (cFacet, pathFraction) = pmeshCI.find_facet(pCoords[pDim:2*pDim], dx, pCellIndex)
                        if cFacet != pmeshCI.NO_FACET:
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
                            mFacet = pmeshCI.cell_entity_index_dict['facet'][pCellIndex][cFacet]
# mFacet should never be the same as the last facet crossed: check this
                            if mFacet == mLastFacet:
                                errorMsg = "%s The mesh index of the facet crossed is %d, the same as the last facet crossed. This should not happen since the particle cannot cross the same facet twice in one move!" % (fncName, mFacet)
                                sys.exit(errorMsg)
                            else:
                                mLastFacet = mFacet
                            # ...and get the value of the facet marker.
# mFacet is a numpy.uint32, but the FacetFunction wants an int argument.
#                            print "type is:", type(mFacet)
                            facValue = pmeshCI.particle_boundary_marker[mFacet]
                            # Check if this facet has a non-zero marker
                            if facValue != 0:
                                # Convert the int marker to its string value...
                                bFlag = str(facValue)
#                                print "bFlag is", bFlag
                                # ...and call the function associated with this value.

# Maybe this needs the whole particle object?  It could be the end-of-the line for the particle.
# Need to check what needs to be in the following arg list:
                                self.pmesh_bcCI.bcFunctionDict[bFlag][sn](psegOut[ipOut], sn, mFacet)

# Not sure if this else is needed:                            else:

#                            # Update the cell index to the new cell
#                                pCellIndex = pmeshCI.cell_neighbor_dict[pCellIndex][facet]
#                                psegOut[ipOut]['cell_index'] = pCellIndex
#                                print "cell_index updated to", psegOut[ipOut]['cell_index']
                            # Update the cell index to the new cell
                            pCellIndex = pmeshCI.cell_neighbor_dict[pCellIndex][cFacet]
                            psegOut[ipOut]['cell_index'] = pCellIndex
                            # If the particle has left the grid, exit this 'while' loop.
                            if pCellIndex == pmeshCI.NO_CELL: break # Exit the 'while' loop
                                
                        else:
                            errorMsg = "%s The cell index of the facet crossed is %d. This should not happen since the particle has left its initial cell cell!" % (fncName, cFacet)
                            sys.exit(errorMsg)
#                       END:if cFacet != pmeshCI.NO_FACET:
#                   END:while not pmeshCI.is_inside(psegOut[ipOut], pCellIndex)

                    # Don't need this since we just look up the cell
                    # index when computing negE above
                    # pseg[ip]['cell_index'] = compute_cell_index(Point(p))

                    # Check that this particle has not been deleted before
                    # incrementing the counter
                    if psegOut[ipOut]['bitflags'] & self.DELETE_FLAG == 0:
                        ipOut += 1
                    else:
                        print "Particle", ipIn, "of species", sn, "has been deleted."

                    # Check if we've reached the end of this segment.  If
                    # so, we need to start writing on a new segment.
                    if (ipOut == self.SEGMENT_LENGTH):
                        particleCount += self.SEGMENT_LENGTH
                        ipOut = 0 # This will cause get_next_out_segment() to be called
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
#        (npSeg, psegIn, psegOut) = psaCI.init_inout_loop()
        (npSeg, psegIn) = psaCI.init_inout_loop()

#        print "move_particles_in_uniform_fields: npSeg, psegIn, psegOut:", npSeg, psegIn, psegOut
        # Get the first segment of particles to move
#        pseg = psaCI.get_next_segment()

        # Loop on segments until a None value is returned
        particleCount = 0
        # ipOut counts through the 'out' array
        ipOut = 0
        while isinstance(psegIn, np_M.ndarray):
#            print "type pseg = ", type(pseg)
#            print pseg['z']
#            print pseg
            # Accelerate this block of particles

            # Loop on the particles in this segment.
            # If a particle has been deleted, skip it.

            # ipIn counts through the 'in' segment
            for ipIn in xrange(npSeg):

                # Skip deleted particles
                if psegIn[ipIn]['bitflags'] & self.DELETE_FLAG != 0: continue

                # Only gets the next 'out' segment if there are more 'in' particles to do
                if ipOut == 0: psegOut = psaCI.get_next_out_segment()

                psegOut[ipOut] = psegIn[ipIn] # This copies everything over, to ensure weights, flags, etc. get copied.

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
                    ipOut = 0 # This will cause get_next_out_segment() to be called
                               # above if there are more 'in' particles to be processed.

    #            print pseg['z']
            (npSeg, psegIn) = psaCI.get_next_segment('in')

        # Set values that will be used to keep track of the particle arrays
#        print "Species", species_name, "npSeg=", npSeg, "psegIn", psegIn, "particleCount", particleCount
        particleCount += ipOut
        psaCI.set_number_of_items('out', particleCount)

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

class ParticleMeshBoundaryConditions_C(object):
    """ParticleMeshBoundaryConditions_C sets up a 2D dictionary
       (bcFunctionDict, indexed by boundary and species) to treat
       kinetic particles incident on a mesh boundary.

       The functions themselves are provided by the user in a
       UserParticleMeshFunctions_C object.
    """

    # Static class variables

    # Particle boundary-conditions are labeled by non-zero bits:
    ABSORB  = 0b1
    REFLECT = 0b1 << 1
    NUMBER_OF_STANDARD_BCs = 2

#    ISEE    = 0b1 << 2        # Ion-stimulated electron emission
#    SEE     = 0b1 << 3        # Secondary-electron emission


# Look for specific boundary conditions

    def __init__(self, speciesNames, pmeshCI, userParticleMeshFunctionsClass, printFlag = False):
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
        # Set local variables from passed parameters

        pBDict = pmeshCI.particle_boundary_dict
        particleBoundaryFlags = pBDict.keys()

        # Initialize a dictionary for the boundary function for each
        # (boundary, species) pair.
        # The first index is the string value of the integer that
        # marks the boundary facets.
        self.bcFunctionDict = dict((bFlag, dict((sp, None) for sp in speciesNames)) for bFlag in particleBoundaryFlags)

        # Get the global default boundary function, if there is one.
        bcFunctionName = 'default_bc'
        if hasattr(userParticleMeshFunctionsClass, bcFunctionName):
            bcGlobalDefaultFunction = getattr(userParticleMeshFunctionsClass, bcFunctionName)
        else:
            bcGlobalDefaultFunction = None

        # Loop on particle boundary names and on particle species
        # names to find the most specific BC.

        # Overwrite the global default function with a function specific to
        # each boundary, if there is one.
        for bFlag in particleBoundaryFlags:
            bcFunctionName = 'default_bc_at_' + pBDict[bFlag]
            if hasattr(userParticleMeshFunctionsClass, bcFunctionName):
                bcBoundaryDefaultFunction = getattr(userParticleMeshFunctionsClass, bcFunctionName)
            else:
                bcBoundaryDefaultFunction = bcGlobalDefaultFunction
            for sp in speciesNames:
                # Overwrite the default for this boundary with a
                # function specific to this boundary and species, if
                # there is one.
                bcFunctionName = 'bc_at_' + pBDict[bFlag] + '_for_' + sp
                if hasattr(userParticleMeshFunctionsClass, bcFunctionName):
                    bcFunction = getattr(userParticleMeshFunctionsClass, bcFunctionName)
                else:
                    bcFunction = bcBoundaryDefaultFunction
                if bcFunction is None:
                    print "ParticleMeshBoundaryConditions_C: No boundary condition specified for", pBDict[bFlag], "/", sp
                elif printFlag:
                    print "ParticleMeshBoundaryConditions_C: Boundary condition for", pBDict[bFlag], "/", sp, "is", bcFunction
                self.bcFunctionDict[bFlag][sp] = bcFunction

        return
#    def __init__(self, particleInputCI, particleCI, printFlag = False):ENDDEF

    def absorb(self):
        return

#    def absorb(self):ENDDEF

    def reflect(self):
        return

#    def reflect(self):ENDDEF

#class ParticleMeshBoundaryConditions_C(object):ENDCLASS
