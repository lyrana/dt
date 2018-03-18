# Particle trajectories

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
__all__ = ['HistoryInput_C',
           'History_C',
#           'FieldHistoryInput_C',
#           'FieldHistory_C',
           'ParticleHistoryInput_C',
           'ParticleHistory_C',
           'TrajectoryInput_C', 
           'Trajectory_C', 
          ]

import sys
import numpy as np_m

# !!! Direct invocation of dolfin for trajectory plots!!!
import dolfin as df_m
import matplotlib.pyplot as mplot_m

#STARTCLASS
class HistoryInput_C(object):
    """An input class for time histories of particle, field, and numerical-parameter
       quantities.

    """

    def __init__(self):
        """Initialize storage for time-histories.

           :cvar max_points: The maximum number of points to record. If 'None', then
                              every point is recorded.

           :cvar scalar_histories: A tuple containing the names of scalar-valued
                                   histories.

           :cvar vector_histories: A tuple containing the names of vector-valued
                                   histories.

        """

        self.max_points = None
        self.scalar_histories = []
        self.vector_histories = []

        return
    
#class HistoryInput_C(object):ENDCLASS

#STARTCLASS
class ParticleHistoryInput_C(object):
    """An input class for time histories of particle-based quantities.

    """

    def __init__(self):
        """Initialize storage for time-histories.

           :cvar max_points: The maximum number of points to record. If 'None', then
                              every point is recorded.

           :cvar scalar_histories: A tuple containing the names of scalar-valued
                                   particle histories.

        """

        self.max_points = None
        self.scalar_histories = []

        return
    
#class ParticleHistoryInput_C(object):ENDCLASS

#STARTCLASS
class FieldHistoryInput_C(object):
    """An input class for time histories of field quantities.

    """

    def __init__(self):
        """Initialize storage for time-histories.

           :cvar max_points: The maximum number of points to record. If 'None', then
                              every point is recorded.

           :cvar scalar_histories: A tuple containing the names of scalar-valued
                                   field histories.

           :cvar vector_histories: A tuple containing the names of vector-valued
                                  field histories.

        """

        self.max_points = None
        self.scalar_histories = []
        self.vector_histories = []

        return
    
#class FieldHistoryInput_C(object):ENDCLASS

#STARTCLASS
class History_C(object):
    """Contains time histories of particle, field, and numerical quantities.

    """

    def __init__(self, histin, ctrl):
        
        """Initialize storage for time-histories.

           :param histin: A HistoryInput_C object.
           :param ctrl: A DTcontrol_C object.

        """

        self.scalar_histories = histin.scalar_histories
        self.ctrl = ctrl

        # Compute a skip interval to stay within the maximum number of
        # datapoints (if specified).
        if histin.max_points is None:
            self.skip = 1
        else:
            self.skip = ctrl.n_timesteps/max_points + 1

        # Length of history-data arrays
        self.npoints = ctrl.n_timesteps/self.skip + 1

        # Initial a counter
        self.counter = 0

        # Create a dictionary to specify the dtype of the history data
        self.history_list_base = ['step', 't'] # These are always recorded
        self.format_list_base = [int, np_m.float32] # Start off the format list with types for
                                             # 'step' and 't'
        
        return
#    def __init__(self, histin, ctrl):ENDDEF


#class History_C(object):
    def record_data(self):
        """Record the history values at this timestep.

        """
        
        
        return    
#    def record_data(self):ENDDEF


#class History_C(object):ENDCLASS


class ParticleHistory_C(History_C):
    """Contains time histories of particle quantities.

    """

    def __init__(self, histin, ctrl, function_dict, species_names):
        
        """Create storage for particle time-histories.

           :param histin: A ParticleHistoryInput_C object.
           :param ctrl: A DTcontrol_C object.
           :param function_dict: The base functions needed to get history data
           :param species_names: A list of particle species names.

           The species names are used to create function calls for the individual
           species from the base functions.

        """

        # Call the parent constructor to complete setting class variables.
        super(self.__class__, self).__init__(histin, ctrl)

        # Particle global histories: for each species, record the history for the sum
        # of all the species, and for each individual species.

        # Add the individual species histories to the user-provided list:
        particleSpeciesHistories = []
#        particleSpeciesHistories[:] = histin.scalar_histories[:]
#        print 'particleSpeciesHistories = ', particleSpeciesHistories
        for hist in histin.scalar_histories:
            for sp in species_names:
                sphist = sp+'_'+hist
                particleSpeciesHistories.append(sphist)
#                species_function_name = 'species_'+hist
#                function_dict[sphist] = (function_dict[species_function_name], sp)
                
        historyList = self.history_list_base + histin.scalar_histories + particleSpeciesHistories
        formatList = self.format_list_base + [np_m.float32 for i in range(len(histin.scalar_histories))] + [np_m.float32 for i in range(len(particleSpeciesHistories))]
#        print 'formatList', formatList

        self.species_scalar_histories = particleSpeciesHistories
        
        # Create the particle-history storage arrays
        self.data_array = np_m.zeros(self.npoints, dtype={'names': historyList, 'formats': formatList})
        
#        # Create the function calls for the individual species


        
        return
#    def __init__(self, histin, ctrl, function_dict, species_names):ENDDEF


#class History_C(object):
    def record_data(self):
        """Record the history values at this timestep.

        """
        
        
        return    
#    def record_data(self):ENDDEF

#class History_C(object):
    def plot(self):
        """Plot the histories owned by this instance of 

        """
        import matplotlib.pyplot as mplot_m

        # Loop on the histories

        tvals = self.data_array['t']
        nlength = self.npoints

        # Plot 'step' vs. 't' for information

        yarray = self.data_array['step']

        plot_title = "timestep number vs. time"
        x_label = "time (s)"
        y_label = "timestep"

        mplot_m.figure()
        mplot_m.plot(tvals, yarray, marker="o")
        mplot_m.title(plot_title)
        mplot_m.xlabel(x_label)
        mplot_m.ylabel(y_label)
        mplot_m.grid(True)
        
        # Plot the requested histories, both total and for each species
        for hist in self.scalar_histories + self.species_scalar_histories:

            yarray = self.data_array[hist]

            plot_title = "%s vs. time" % (hist)
            x_label = "time (s)"
            y_label = hist

            mplot_m.figure()
#            mplot_m.plot(tvals[0:nlength], data_arr[comp][0:nlength], marker="o")
            mplot_m.plot(tvals, yarray, marker="o")
            mplot_m.title(plot_title)
            mplot_m.xlabel(x_label)
            mplot_m.ylabel(y_label)
            mplot_m.grid(True)
            
        mplot_m.show()
            
        return
#    def plot(self):ENDDEF



#class History_C(object):ENDCLASS




#STARTCLASS
class TrajectoryInput_C(object):
    """Parameters that specify particle trajectory
    """

    def __init__(self):
        """ Initialize variables
        """

        # If max_points is specified, no more than this number of points will be saved.
        self.max_points = None

        # The number of points that will be saved. The upper limit is the number of
        # timesteps plus extra_points.
        self.npoints = None

        # An array length equal to the number of timesteps+1 may not be enough to record all
        # points on the trajectory, since boundary-crossings are always recorded.  To allow for
        # one boundary-crossing where the particle is absorbed at the boundary,
        # extra_points can be set to 1.  For a particle reflecting multiple times off a
        # boundary, a larger value can be used.
        self.extra_points = None

        self.explicit_dict = None
        self.implicit_dict = None
        self.neutral_dict = None

        return

#class TrajectoryInput_C(object):ENDCLASS

#STARTCLASS
class Trajectory_C(object):
    """Trajectory_C has a collection of particle trajectories indexed
       by species name.

       Particles that you want a trajectory for are marked with
       TRAJECTORY_FLAG, and are thereafter identifed by their index in
       the particle storage array.

       Storage is created by 'create_trajectory()' when particles with
       the TRAJECTORY_FLAG turned on are encountered.
    """

    ### Static class variables

    # Value signaling that the particle no longer exists.  E.g., the particle hit a
    # wall and was deleted.  Trajectory data for the particle exist up to the point
    # where it was deleted.
    NO_PINDEX = -1

    def __init__(self, trajin, ctrl, explicit_species, implicit_species, neutral_species):
        """Set up the initial trajectory storage.  

           :var trajin.max_points: Maximum number of trajectory data points to be
                                   collected.  If set to 'None', then every step on the
                                   trajectory is saved.

           :var npoints: The number of datapoints that will be saved in a trajectory.  If
                         max_points is set, then npoints = max_points.  If max_points is
                         'None', then npoints = n_timesteps, i.e., every step is saved.

           :var particle_index_list: particle_index_list[sp] is a list of the
                                     particle indices of species sp that were selected for
                                     trajectories.

           :var data_list: data_list[sp][itraj] is a numpy array of length npoints that
                           stores the trajectory data for trajectory itraj of species sp.

           :var explicit_dict: a dictionary in numpy dtype format giving the string names
                               of particle attributes that can be saved in a trajectory.

           The storage uses numpy arrays of length self.npoints           
           npoints is at most the number of timesteps.

           TrajDict is a dictionary of the trajectory variables.
           e.g.,  {'names': ['x', 'y'], 'formats': [np_m.float64, np_m.float64]}

        """

        self.last_step = None # This variable can be used to check if the simulation
                              # step in the current call is the same as the step in
                              # the last call, and avoid storing the same data twice.
                              # If the particle is deleted just after being pushed,
                              # there's special handling to allow its last data to be
                              # recorded.

        # Compute a skip interval to stay within the maximum number of
        # datapoints (if specified).
        if trajin.max_points is None:
            self.skip = 1
        else:
            self.skip = ctrl.n_timesteps/trajin.max_points + 1

        # Length of trajectory data arrays
        # If the location of boundary-crossings is recorded, then you can run out of space
        self.npoints = ctrl.n_timesteps/self.skip + 1 + trajin.extra_points

        # Need these to get the right trajectory variables
        self.explicit_species = explicit_species
        self.implicit_species = implicit_species
        self.neutral_species = neutral_species

        # The dictionary is in the form of a numpy dtype, giving names
        # and types of the trajectory values, as specified by the user's input.
        self.explicit_dict = trajin.explicit_dict
        self.implicit_dict = trajin.implicit_dict
        self.neutral_dict = trajin.neutral_dict

        # For each species, there is
        #    1. a list of trajectory-particle IDs
        #    2. a list of trajectory data arrays, one for each marked particle.
        #    3. a list of trajectory lengths, one for each marked particle.
        # These lists are indexed by species name, so they are dictionaries:
        self.particle_index_list = {} # List of particle storage indices
        self.particle_unique_id_list = {} # List of particle unique_IDs
        self.data_list = {}
        self.trajectory_length = {}

        for sp in self.explicit_species + self.implicit_species + self.neutral_species:
            self.particle_index_list[sp] = []
            self.particle_unique_id_list[sp] = []
            self.data_list[sp] = []
            self.trajectory_length[sp] = []

        return
#    def __init__(self, trajin, ctrl, explicit_species, implicit_species, neutral_species):ENDDEF

    def create_trajectory(self, species_name, pindex, dynamics_type, unique_id_int=None):
        """Add storage for another particle trajectory.

           :param species_name: Name of the species that the marked
                                particle belongs to.
           :param pindex: Index of the marked particle in the particle
                          storage array.
           :param int unique_id_int: The unique ID integer of this particle.

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        # Add the particle's index to the list of marked particles for this species.
        self.particle_index_list[species_name].append(pindex)
        self.particle_unique_id_list[species_name].append(unique_id_int)

#        print "create_traj: The current particle indices are:", self.particle_index_list[species_name]
        
        self.trajectory_length[species_name].append(0) # Set the trajectory point-count for this
                                                      # particle to 0

        # Add a numpy array for this particles's trajectory data
        if dynamics_type == 'explicit':
            self.data_list[species_name].append(np_m.empty(self.npoints, dtype=self.explicit_dict))
        elif dynamics_type == 'implicit':
            self.data_list[species_name].append(np_m.empty(self.npoints, dtype=self.implicit_dict))
        elif dynamics_type == 'neutral':
            self.data_list[species_name].append(np_m.empty(self.npoints, dtype=self.neutral_dict))
        else:
            errorMsg = "(DnT ERROR) %s\tdynamics_type %s is unknown for species %s" % (fncName, dynamics_type, species_name)
            raise RuntimeError(errorMsg)

        return
#    def create_trajectory(self, species_name, dynamics_type):ENDDEF

#class Trajectory_C(object):
    def plot(self):
        """Plot the accumulated particle trajectories.

        """
        import matplotlib.pyplot as mplot_m

#        tmin = 0.0
#        tmax = tmin + (self.npoints-1)*self.data_interval
#        tvals = np_m.linspace(tmin, tmax, self.npoints)


#        print 'tval = ', tvals

        # Temporary array for plotting
#        yvals = np_m.empty(self.npoints)

        for sp in self.explicit_species + self.implicit_species + self.neutral_species:
            if len(self.data_list[sp]) == 0: continue
            comps = self.data_list[sp][0][0].dtype.names[2:] # Skip the 'step' and 't' fields
#            print 'comps =', comps
            for it in xrange(len(self.particle_index_list[sp])):
                if self.particle_unique_id_list[sp][it] is None:
                    ip = self.particle_index_list[sp][it]
                    plot_title = "%s: Traj# %d Particle id %d" % (sp, it, ip)
                else:
                    ip = self.particle_unique_id_list[sp][it]
                    plot_title = "%s: Traj# %d Particle uid %d" % (sp, it, ip)
                data_arr = self.data_list[sp][it]
                tvals = self.data_list[sp][it]['t']
                nlength = self.trajectory_length[sp][it]
#                print 'tvals = ', tvals[0:nlength]

                # Plots vs. time
                
                x_label = "time (s)"
                for comp in comps:
                    mplot_m.figure()
                    mplot_m.plot(tvals[0:nlength], data_arr[comp][0:nlength], marker="o")
                    mplot_m.title(plot_title)
                    mplot_m.xlabel(x_label)
                    mplot_m.ylabel(comp)
                    mplot_m.grid(True)
#                    mplot_m.show()

                # Phase-space plots
                
                # x vs. y
                if 'x' in comps:
                    if 'y' in comps:
                        mplot_m.figure()
                        mplot_m.plot(data_arr['x'][0:nlength], data_arr['y'][0:nlength], marker="o")
                        mplot_m.title(plot_title)
                        mplot_m.xlabel('x')
                        mplot_m.ylabel('y')
                        mplot_m.grid(True)
#                        mplot_m.show()
                    if 'ux' in comps:
                        mplot_m.figure()
                        mplot_m.plot(data_arr['x'][0:nlength], data_arr['ux'][0:nlength], marker="o")
                        mplot_m.title(plot_title)
                        mplot_m.xlabel('x')
                        mplot_m.ylabel('ux')
                        mplot_m.grid(True)
#                        mplot_m.show()
                
                mplot_m.show()
        return
#    def plot(self):ENDDEF

#class Trajectory_C(object):
    def plot_trajectories_on_mesh(self, mesh, plot_title, hold_plot=False):
        """Plot the (x, y) trajectory points on top of a mesh plot
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

        # Create a plotter to display the trajectory.
        # !!! Direct call to dolfin!!!
        plotter=df_m.plot(mesh, title=plot_title)

        for sp in self.explicit_species + self.implicit_species + self.neutral_species:
            if len(self.data_list[sp]) == 0: continue
            comps = self.data_list[sp][0][0].dtype.names[2:] # Skip the 'step' and 't' fields
            # Loop on trajectories for this species
            for it in xrange(len(self.particle_index_list[sp])):
                nlength = self.trajectory_length[sp][it]
                if nlength == 1:
                    print fncName, "*** DT Warning: Trajectory", it, "for species", sp, "has only one point. Not plotting it.***"
                    continue
                data_arr = self.data_list[sp][it]
#                print "data_arr[x] =", data_arr['x'][0:nlength]
                # x vs. y
                if 'x' in comps:
                    if 'y' in comps:
                        mplot_m.plot(data_arr['x'][0:nlength], data_arr['y'][0:nlength], marker="o")
# VTK plotter:                        
#                        path = np_m.empty(2*nlength, dtype=np_m.float64) # NB: dtype has to be double for add_polygon()
#                        path[0::2] = data_arr['x'][0:nlength] # Start at 0, with stride 2
#                        path[1::2] = data_arr['y'][0:nlength] # Start at 1, with stride 2
#                        print fncName, 'path =', path
#                        print fncName, 'Trajectory', it, 'has ', path.size, 'points'
                        # Replace with a point plot:
                        # plotter.add_polygon(path)

        mplot_m.show() # Display the plot

        return
#    def plot_trajectories_on_mesh(self):ENDDEF

#class Trajectory_C(object):ENDCLASS

