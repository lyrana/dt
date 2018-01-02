# Particle trajectories

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
__all__ = ['TrajectoryInput_C', 
           'Trajectory_C', 
          ]

import sys
import numpy as np_M

# !!! Direct invocation of dolfin for trajectory plots!!!
import dolfin as df_M

#STARTCLASS
class TrajectoryInput_C(object):
    """Parameters that specify particle trajectory
    """

    def __init__(self):
        """ Initialize variables
        """

        # Upper limit is the number of timesteps
        self.maxpoints = None

        self.npoints = None

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

    def __init__(self, trajinCI, ctrl, explicit_species, implicit_species, neutral_species):
        """Set up the initial storage.  

           :var trajinCI.maxpoints: Maximum number of trajectory data points to be
                                    collected.  If set to 'None', then every step on the
                                    trajectory is saved.

           :var npoints: The number of datapoints that will be saved in a trajectory.  If
                         maxpoints is set, then npoints = maxpoints.  If maxpoints is
                         'None', then npoints = n_timesteps, i.e., every step is saved.

           :var ParticleIdList: ParticleIdList[sp] is a list of the particle indices of
                                species sp that were selected for trajectories.

           :var DataList: DataList[sp][itraj] is a numpy array of length npoints that
                          stores the trajectory data for trajectory itraj of species sp.

           :var explicit_dict: a dictionary in numpy dtype format giving the string names
                              of particle attributes that can be saved in a trajectory.

           The storage uses numpy arrays of length self.npoints           
           npoints is at most the number of timesteps.

           TrajDict is a dictionary of the trajectory variables.
           e.g.,  {'names': ['x', 'y'], 'formats': [np_M.float64, np_M.float64]}
        """

        # Counter for the number of datapoints stored
        self.count = 0 # XX Not used?

        self.last_step = None # This variable can be used to check if the simulation
                              # step in the current call is the same as the step in
                              # the last call, and avoid storing the same data twice.
                              # If the particle is deleted just after being pushed,
                              # there's special handling to allow its last data to be
                              # recorded.

        # Compute a skip interval to stay within the maximum number of
        # datapoints (if specified).
        if trajinCI.maxpoints == None:
            self.skip = 1
        else:
            self.skip = ctrl.n_timesteps/trajinCI.maxpoints + 1

        # Length of trajectory data arrays
#        self.npoints = ctrl.n_timesteps/self.skip + 2 # Add 2 because a newly-created trajectory particle automatically gets recorded, and may then be recorded again 
        self.npoints = ctrl.n_timesteps/self.skip + 1

        # Time interval between datapoints
        self.data_interval = self.skip * ctrl.dt

        # Need these to get the right trajectory variables
        self.explicit_species = explicit_species
        self.implicit_species = implicit_species
        self.neutral_species = neutral_species

        # The dictionary is in the form of a numpy dtype, giving names
        # and types of the trajectory values, as specified by the user's input.
        self.explicit_dict = trajinCI.explicit_dict
        self.implicit_dict = trajinCI.implicit_dict
        self.neutral_dict = trajinCI.neutral_dict

        # For each species, there is
        #    1. a list of trajectory-particle IDs
        #    2. a list of trajectory data arrays, one for each marked particle.
        #    3. a list of trajectory lengths, one for each marked particle.
        # These lists are indexed by species name, so they are dictionaries:
        self.ParticleIdList = {} # List of particle IDs
        self.DataList = {}
        self.TrajectoryLength = {}

        for sp in self.explicit_species + self.implicit_species + self.neutral_species:
            self.ParticleIdList[sp] = []
            self.DataList[sp] = []
            self.TrajectoryLength[sp] = []

        return
#    def __init__(self, trajinCI, ctrl, explicit_species, implicit_species, neutral_species):ENDDEF

    def create_trajectory(self, species_name, pindex, dynamics_type):
        """Add storage for another particle trajectory.

           :param species_name: Name of the species that the marked
                                particle belongs to.
           :param pindex: Index of the marked particle in the particle
                          storage array.


        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        self.ParticleIdList[species_name].append(pindex)

#        print "create_traj: The current particle indices are:", self.ParticleIdList[species_name]
        
        self.TrajectoryLength[species_name].append(0) # Set the count to 0

        # Add a numpy array for the trajectory data
        if dynamics_type == 'explicit':
            self.DataList[species_name].append(np_M.empty(self.npoints, dtype=self.explicit_dict))
        elif dynamics_type == 'implicit':
            self.DataList[species_name].append(np_M.empty(self.npoints, dtype=self.implicit_dict))
        elif dynamics_type == 'neutral':
            self.DataList[species_name].append(np_M.empty(self.npoints, dtype=self.neutral_dict))
        else:
            errorMsg = "(DnT ERROR) %s\tdynamics_type %s is unknown for species %s" % (fncName, dynamics_type, species_name)
            raise RuntimeError(errorMsg)

        return
#    def create_trajectory(self, species_name, dynamics_type):ENDDEF

#class Trajectory_C(object):
    def plot_trajectories(self):
        """Plot the nth trajectory
        """
        import matplotlib.pyplot as mplot_m

#        tmin = 0.0
#        tmax = tmin + (self.npoints-1)*self.data_interval
#        tvals = np_M.linspace(tmin, tmax, self.npoints)


#        print 'tval = ', tvals

        # Temporary array for plotting
#        yvals = np_M.empty(self.npoints)

        for sp in self.explicit_species + self.implicit_species + self.neutral_species:
#            self.ParticleIdList[sp] = []
#            self.DataList[sp] = []
            comps = self.DataList[sp][0][0].dtype.names[2:] # Skip the 'step' and 't' fields
#            print 'comps =', comps
            for it in xrange(len(self.ParticleIdList[sp])):
                ip = self.ParticleIdList[sp][it]
                data_arr = self.DataList[sp][it]
                tvals = self.DataList[sp][it]['t']
                nlength = self.TrajectoryLength[sp][it]
#                print 'tvals = ', tvals[0:nlength]
                # Plots vs. time
                for comp in comps:
                    mplot_m.plot(tvals[0:nlength], data_arr[comp][0:nlength])
                    mplot_m.title("%s: Traj# %d Particle id %d" % (sp, it, ip))
                    mplot_m.xlabel('time (s)')
                    mplot_m.ylabel(comp)
                    mplot_m.grid(True)
                    mplot_m.show()
                # Phase-space plots
                # x vs. y
                if 'x' in comps:
                    if 'y' in comps:
                        mplot_m.plot(data_arr['x'][0:nlength], data_arr['y'][0:nlength])
                        mplot_m.xlabel('x')
                        mplot_m.ylabel('y')
                        mplot_m.grid(True)
                        mplot_m.show()
                    if 'ux' in comps:
                        mplot_m.plot(data_arr['x'][0:nlength], data_arr['ux'][0:nlength])
                        mplot_m.xlabel('x')
                        mplot_m.ylabel('ux')
                        mplot_m.grid(True)
                        mplot_m.show()
                
        return
#    def plot_trajectories(self):ENDDEF

#class Trajectory_C(object):
    def plot_trajectories_on_mesh(self, mesh, plot_title, hold_plot=False):
        """Plot the (x, y) trajectory points on top of a mesh plot
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

        # Create a plotter to display the trajectory.
        # !!! Direct call to dolfin!!!
        plotter=df_M.plot(mesh, title=plot_title)

        for sp in self.explicit_species + self.implicit_species + self.neutral_species:
            comps = self.DataList[sp][0][0].dtype.names[2:] # Skip the 'step' and 't' fields
            # Loop on trajectories for this species
            for it in xrange(len(self.ParticleIdList[sp])):
#                ip = self.ParticleIdList[sp][it] # Not used
                nlength = self.TrajectoryLength[sp][it]
                if nlength == 1:
                    print fncName, "*** DT Warning: Trajectory", it, "for species", sp, "has only one point. Not plotting it.***"
                    continue
                data_arr = self.DataList[sp][it]
#                print "data_arr[x] =", data_arr['x'][0:nlength]
                # x vs. y
                if 'x' in comps:
                    if 'y' in comps:
                        path = np_M.empty(2*nlength, dtype=np_M.float64) # NB: dtype has to be double for add_polygon()
                        path[0::2] = data_arr['x'][0:nlength] # Start at 0, with stride 2
                        path[1::2] = data_arr['y'][0:nlength] # Start at 1, with stride 2
#                        print fncName, 'path =', path
#                        print fncName, 'Trajectory', it, 'has ', path.size, 'points'
                        plotter.add_polygon(path)

        plotter.plot()
        if hold_plot is True: df_M.interactive() # Stops the plot from disappearing

        return
#    def plot_trajectories_on_mesh(self):ENDDEF

#class Trajectory_C(object):ENDCLASS
