# Particle trajectories

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
__all__ = ['', 
           '', 
           '', ]

import sys
import numpy as np_M

# !!! Direct invocation of dolfin for trajectory plots!!!
import dolfin as df_M

class Trajectory_C(object):
    """Trajectory_C has a collection of particle trajectories for one species.
    """

    def __init__(self, trajinCI, ctrlCI, explicit_species, implicit_species, neutral_species):
        """Set up the initial storage.  
           TrajLength is at most the number of timesteps.
           TrajDict is a dictionary of the trajectory variables.
           e.g.,  {'names': ['x', 'y'], 'formats': [np_M.float64, np_M.float64]}
        """

        # Counter for the number of datapoints stored
        self.count = 0

        # Compute a skip interval to stay within the maximum number of
        # datapoints (if specified).
        if trajinCI.maxpoints == None:
            self.skip = 1
        else:
            self.skip = ctrlCI.nsteps/trajinCI.maxpoints + 1

        # Length of trajectory data arrays
        self.npoints = ctrlCI.nsteps/self.skip + 1

        # Time interval between datapoints
        self.data_interval = self.skip * ctrlCI.dt

        # Need these to get the right trajectory variables
        self.explicit_species = explicit_species
        self.implicit_species = implicit_species
        self.neutral_species = neutral_species

        # The dictionary is in the form of a numpy dtype, giving names
        # and types of the trajectory values, as specified by the user's input.
        self.explicitDict = trajinCI.explicitDict
        self.implicitDict = trajinCI.implicitDict
        self.neutralDict = trajinCI.neutralDict

        # For each species, there is
        #    1. a list of trajectory-particle IDs
        #    2. a list of trajectory data arrays, one for each marked particle.
        # These lists are indexed by species name, so they are dictionaries:
        self.ParticleIdList = {}
        self.DataList = {}

        for sp in self.explicit_species + self.implicit_species:
            self.ParticleIdList[sp] = []
            self.DataList[sp] = []

        # Create a ndarray for 1 particle that we can use like a SA.
#            self.explicit_1particle_arr = np_M.empty(1, dtype=self.explicitDict))
#            self.implicit_1particle_arr = np_M.empty(1, dtype=self.implicitDict))

        return

    def create_trajectory(self, species_name, dynamics_type):
        """Add storage for another particle trajectory.
        """

        # Add a numpy array for the data
        if dynamics_type == 'explicit':
            self.DataList[species_name].append(np_M.empty(self.npoints, dtype=self.explicitDict))
        elif dynamics_type == 'implicit':
            self.DataList[species_name].append(np_M.empty(self.npoints, dtype=self.implicitDict))
        elif dynamics_type == 'neutral':
            self.DataList[species_name].append(np_M.empty(self.npoints, dtype=self.neutralDict))
        else:
            error_msg = "Trajectory_C:create_trajectory: dynamics_type %s is unknown" % dynamics_type
            sys.exit(error_msg)

        return

#class Trajectory_C(object):
    def plot_trajectories(self):
        """Plot the nth trajectory
        """
        import matplotlib.pyplot as mplot_M

        tmin = 0.0
        tmax = tmin + (self.npoints-1)*self.data_interval
        tvals = np_M.linspace(tmin, tmax, self.npoints)
#        print 'tval = ', tvals

        # Temporary array for plotting
#        yvals = np_M.empty(self.npoints)

        for sp in self.explicit_species + self.implicit_species + self.neutral_species:
#            self.ParticleIdList[sp] = []
#            self.DataList[sp] = []
            comps = self.DataList[sp][0][0].dtype.names
#            print 'comps =', comps
            for it in xrange(len(self.ParticleIdList[sp])):
                ip = self.ParticleIdList[sp][it]
                data_arr = self.DataList[sp][it]
                # Plots vs. time
                for comp in comps:
                    mplot_M.plot(tvals, data_arr[comp])
                    mplot_M.title("%s: Traj# %d Particle id %d" % (sp, it, ip))
                    mplot_M.xlabel('time (s)')
                    mplot_M.ylabel(comp)
                    mplot_M.grid(True)
                    mplot_M.show()
                # Phase-space plots
                # x vs. y
                if 'x' in comps:
                    if 'y' in comps:
                        mplot_M.plot(data_arr['x'], data_arr['y'])
                        mplot_M.xlabel('x')
                        mplot_M.ylabel('y')
                        mplot_M.grid(True)
                        mplot_M.show()
                    if 'ux' in comps:
                        mplot_M.plot(data_arr['x'], data_arr['ux'])
                        mplot_M.xlabel('x')
                        mplot_M.ylabel('ux')
                        mplot_M.grid(True)
                        mplot_M.show()
                
        return
#    def plot_trajectories(self):ENDDEF

#class Trajectory_C(object):
    def plot_trajectories_on_mesh(self, mesh):
        """Plot the trajectories on top of a mesh plot
        """

        # Create a plotter to display the trajectory.
        # !!! Direct call to dolfin !!!
        plotter=df_M.plot(mesh, title="Trajectory plot")

        for sp in self.explicit_species + self.implicit_species + self.neutral_species:
            comps = self.DataList[sp][0][0].dtype.names
            # Loop on trajectories for this species
            for it in xrange(len(self.ParticleIdList[sp])):
                ip = self.ParticleIdList[sp][it]
                data_arr = self.DataList[sp][it]
                # x vs. y
                if 'x' in comps:
                    if 'y' in comps:
                        path = np_M.empty(2*data_arr['x'].size, dtype=np_M.float64) # dtype has to be double for add_polygon()
                        path[0::2] = data_arr['x']
                        path[1::2] = data_arr['y']
#                        print 'path =', path
                        plotter.add_polygon(path)

        plotter.plot()
#        df_M.interactive() # Stops the plot from disappearing
                
        return
#    def plot_trajectories_on_mesh(self):ENDDEF

#class Trajectory_C(object):ENDCLASS
