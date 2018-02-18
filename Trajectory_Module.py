# Particle trajectories

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
__all__ = ['TrajectoryInput_C', 
           'Trajectory_C', 
          ]

import sys
import numpy as np_m

# !!! Direct invocation of dolfin for trajectory plots!!!
import dolfin as df_m
import matplotlib.pyplot as mplot_m

#STARTCLASS
class TrajectoryInput_C(object):
    """Parameters that specify particle trajectory
    """

    def __init__(self):
        """ Initialize variables
        """

        # If maxpoints is specified, no more than this number of points will be saved.
        self.maxpoints = None

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
        """Set up the initial storage.  

           :var trajin.maxpoints: Maximum number of trajectory data points to be
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
        if trajin.maxpoints == None:
            self.skip = 1
        else:
            self.skip = ctrl.n_timesteps/trajin.maxpoints + 1

        # Length of trajectory data arrays
        # If the location of boundary-crossings is recorded, then you can run out of space
        self.npoints = ctrl.n_timesteps/self.skip + 1 + trajin.extra_points

        # Time interval between datapoints
        self.data_interval = self.skip * ctrl.dt

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
        self.ParticleIdList = {} # List of particle storage indices
        self.ParticleUniqueIdList = {} # List of particle unique_IDs
        self.DataList = {}
        self.TrajectoryLength = {}

        for sp in self.explicit_species + self.implicit_species + self.neutral_species:
            self.ParticleIdList[sp] = []
            self.ParticleUniqueIdList[sp] = []
            self.DataList[sp] = []
            self.TrajectoryLength[sp] = []

        return
#    def __init__(self, trajin, ctrl, explicit_species, implicit_species, neutral_species):ENDDEF

    def create_trajectory(self, species_name, pindex, dynamics_type, unique_id_int=None):
        """Add storage for another particle trajectory.

           :param species_name: Name of the species that the marked
                                particle belongs to.
           :param pindex: Index of the marked particle in the particle
                          storage array.
           :param int unique_ID: The unique ID integer of this particle.

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        # Add the particle's index to the list of marked particles for this species.
        self.ParticleIdList[species_name].append(pindex)
        self.ParticleUniqueIdList[species_name].append(unique_id_int)

#        print "create_traj: The current particle indices are:", self.ParticleIdList[species_name]
        
        self.TrajectoryLength[species_name].append(0) # Set the trajectory point-count for this
                                                      # particle to 0

        # Add a numpy array for this particles's trajectory data
        if dynamics_type == 'explicit':
            self.DataList[species_name].append(np_m.empty(self.npoints, dtype=self.explicit_dict))
        elif dynamics_type == 'implicit':
            self.DataList[species_name].append(np_m.empty(self.npoints, dtype=self.implicit_dict))
        elif dynamics_type == 'neutral':
            self.DataList[species_name].append(np_m.empty(self.npoints, dtype=self.neutral_dict))
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
#        tvals = np_m.linspace(tmin, tmax, self.npoints)


#        print 'tval = ', tvals

        # Temporary array for plotting
#        yvals = np_m.empty(self.npoints)

        for sp in self.explicit_species + self.implicit_species + self.neutral_species:
            if len(self.DataList[sp]) == 0: continue
            comps = self.DataList[sp][0][0].dtype.names[2:] # Skip the 'step' and 't' fields
#            print 'comps =', comps
            for it in xrange(len(self.ParticleIdList[sp])):
                if self.ParticleUniqueIdList[sp][it] is None:
                    ip = self.ParticleIdList[sp][it]
                    plot_title = "%s: Traj# %d Particle id %d" % (sp, it, ip)
                else:
                    ip = self.ParticleUniqueIdList[sp][it]
                    plot_title = "%s: Traj# %d Particle uid %d" % (sp, it, ip)
                data_arr = self.DataList[sp][it]
                tvals = self.DataList[sp][it]['t']
                nlength = self.TrajectoryLength[sp][it]
#                print 'tvals = ', tvals[0:nlength]
                # Plots vs. time
                x_label = "'time (s)"
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
#    def plot_trajectories(self):ENDDEF

#class Trajectory_C(object):
    def plot_trajectories_on_mesh(self, mesh, plot_title, hold_plot=False):
        """Plot the (x, y) trajectory points on top of a mesh plot
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

        # Create a plotter to display the trajectory.
        # !!! Direct call to dolfin!!!
        plotter=df_m.plot(mesh, title=plot_title)

        for sp in self.explicit_species + self.implicit_species + self.neutral_species:
            if len(self.DataList[sp]) == 0: continue
            comps = self.DataList[sp][0][0].dtype.names[2:] # Skip the 'step' and 't' fields
            # Loop on trajectories for this species
            for it in xrange(len(self.ParticleIdList[sp])):
                nlength = self.TrajectoryLength[sp][it]
                if nlength == 1:
                    print fncName, "*** DT Warning: Trajectory", it, "for species", sp, "has only one point. Not plotting it.***"
                    continue
                data_arr = self.DataList[sp][it]
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
