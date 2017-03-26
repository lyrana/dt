# Make a mesh that DT can use

#__version__ = 0.1
#__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

"""UserMesh_y_Fields defines the mesh, marks boundaries for field and
   particle boundary-conditions (BCs), and creates the field solvers.
   The type of mesh may be FE or FD.  In this case, it's easiest to
   apply the BCs while the mesh is being created, since at first, the
   boundaries are Cartesian coordinate lines that are easy to specify.
   The mesh is then stretched and transformed to a cylindrical shape.
"""

import sys
import os
import numpy as np_M
import math
# !!! Direct invocation of dolfin. OK because UserMesh_C is a
# sub-class of Mesh_C !!!
import dolfin as df_M

from Dolfin_Module import Mesh_C
from Dolfin_Module import PoissonSolve_C

import UserUnits_Module as U_M

# Create definitions of functions to test if points are on the boundaries
# The XBoundary function takes care of the two Dirichlet boundaries, and
# set_all() marks every facet with a 2.

# Don't actually need the AllFacets function for this problem
class AllFacets(df_M.SubDomain):
    """The AllFacets class is a specialized SubDomain
    """
    # Return true for every location x that is on a boundary
    def inside(self, x, on_boundary):
        return on_boundary

class XBoundary(df_M.SubDomain):
    """The XBoundary class is a specialized SubDomain
    """
    # Set the X value of the boundary in the constructor
    def __init__(self, x_boundary_value, map_tol=1.0e-10):
        self.x_boundary_value=x_boundary_value
        # Pass the map_tol argument on to the base class
        super(self.__class__, self).__init__(map_tol=1.0e-10)

    # Return true for points inside the subdomain
    def inside(self, x, on_boundary):
        tol = 1.0e-10
        return on_boundary and abs(x[0]-self.x_boundary_value) < tol

class YBoundary(df_M.SubDomain):
    """The YBoundary class is a specialized SubDomain
    """
    # Set the X value of the boundary in the constructor
    def __init__(self, y_boundary_value, map_tol=1.0e-10):
        self.y_boundary_value=y_boundary_value
        # Pass the map_tol argument on to the base class
        super(self.__class__, self).__init__(map_tol=1.0e-10)

    # Return true for points inside the subdomain
    def inside(self, x, on_boundary):
        """x[0|1|2] contains the x, y, z coordinates of a point.
        """
        tol = 1.0e-10
        return on_boundary and abs(x[1]-self.y_boundary_value) < tol

# User exposes whatever mesh parameters are useful in __init__ and
# these can be set in __main__

class UserMesh_C(Mesh_C):
    """UserMesh_C is derived from Mesh_C, and includes marked
       subdomains (e.g., Dirichlet boundaries).  It is to be edited by
       the user to specify the simulation mesh.

       The units are MKS by default (i.e., if no conversion
       factor is applied), but any units can be used provided the
       conversion to MKS is available in the UserUnits_M.py module.
    """

    # Select the unit system to be used for input parameters.
    Convert = U_M.MyPlasmaUnits_C

    def __init__(self, meshInputCI=None, compute_dictionaries=False, compute_tree=False, plot_flag=False):
        """
            The class UserMesh_C contains these attributes:
                1. A mesh.
                2. A fieldBoundaryMarker (MeshFunction) that marks
                   mesh boundaries for field boundary-conditions.
                3. A particleBoundaryMarker (MeshFunction) that marks
                   mesh boundaries for particle boundary-conditions.
        """

        if meshInputCI.mesh_file is None:
            self.create_mesh(meshInputCI, plot_flag)
            # don't need another mesh plot
            plot_flag = False

        # Call the parent constructor to complete setting class variables.
        mesh_file = meshInputCI.mesh_file
        super(self.__class__, self).__init__(mesh_file=mesh_file, compute_dictionaries=compute_dictionaries, compute_tree=compute_tree, plot_flag=plot_flag)
        
        self.field_boundary_dict = meshInputCI.field_boundary_dict
        self.particle_boundary_dict = meshInputCI.particle_boundary_dict

        # Read boundary markers from file, if provided
        fieldBoundaryFile = meshInputCI.field_boundary_file
        particleBoundaryFile = meshInputCI.particle_boundary_file
        if meshInputCI.field_boundary_file is not None:
            fieldBoundaryMarker = df_M.MeshFunction("size_t", self.mesh, fieldBoundaryFile)
            self.field_boundary_marker = fieldBoundaryMarker
        if particleBoundaryFile is not None:
            particleBoundaryMarker = df_M.MeshFunction("size_t", self.mesh, particleBoundaryFile)
# or equivalently:       particleBoundaryMarker = df_M.MeshFunctionSizet(self.mesh, "Pbcs_quarter_circle_mesh_crossed.xml")
            self.particle_boundary_marker = particleBoundaryMarker

        return  
#    def __init__(self, meshInputCI=None, compute_dictionaries=False, compute_tree=False, plot_flag=False):ENDDEF

# Inherited from Mesh_C:
#    def copy(self):
#        return copy.deepcopy(self)

#class UserMesh_C(Mesh_C):
    def create_mesh(self, meshInputCI, plot_flag):
        """Create a mesh and mark subdomains (e.g., Dirichlet
           boundaries and particle boundaries) according to the user's
           specifications.
        """

        # Radial limits of mesh
        rmin = meshInputCI.rmin
        rmax = meshInputCI.rmax

        stretch = meshInputCI.stretch
        nr = meshInputCI.nr

        # The final mesh will go from 0 to tmax:
        tmax = meshInputCI.tmax
        nt = meshInputCI.nt

        # Options: 'left, 'right', 'left/right', 'crossed'
        diagonal = meshInputCI.diagonal

        # Boundary conditions for fields and particles
        fieldBoundaryDict = meshInputCI.field_boundary_dict
        particleBoundaryDict = meshInputCI.particle_boundary_dict

# First, make a rectangular mesh

#       args: RectangleMesh(x0, y0, x1, y1, nx, ny, diagonal="right")

# Switch these for CEE/Laptop
        if df_M.DOLFIN_VERSION_STRING > "1.5.0":
# v > 1.5:
            mesh = df_M.RectangleMesh(df_M.Point(rmin, 0.0), df_M.Point(rmax, 1.0), nr, nt, diagonal)
        else:
# v = 1.5:
            mesh = df_M.RectangleMesh(rmin, 0.0, rmax, 1.0, nr, nt, diagonal)

#V = df.FunctionSpace(mesh, 'Lagrange', 1)

        # Field boundary conditions

        if fieldBoundaryDict is not None:
            # Set the Dirichlet boundary-conditions on the Rectangle,
            # before transforming it, since it's a simple rectangle then.

            # Create a MeshFunction object that contains the mesh, and is
            # defined on mesh facets.  The function has a (size_t) value
            # of 0, 1, or 2, to identify the boundary-condition applied.
            # The function values are set later.

            # A boundary has dimension 1 less than the domain:
    # equivalent:       fieldBoundaryMarker = df_M.MeshFunction('size_t', mesh, mesh.topology().dim()-1)
            fieldBoundaryMarker = df_M.FacetFunction('size_t', mesh)

            # Initialize all mesh facets with a default value of 0
            # Then overwrite boundary facets that have Dirichlet BCs.
            fieldBoundaryMarker.set_all(0)

            # Create mesh subdomains to apply boundary-conditions
            Gamma_rmin = XBoundary(rmin) # Gamma_rmin is the lower radial boundary
            Gamma_rmax = XBoundary(rmax) # Gamma_rmax is the upper radial boundary

            # Mark these subdomains (boundaries) with non-zero integers
            rmin_indx = fieldBoundaryDict['rmin']
            rmax_indx = fieldBoundaryDict['rmax']
            Gamma_rmin.mark(fieldBoundaryMarker, rmin_indx) # Mark the lower radial
                                                       # boundary with rmin_indx
            Gamma_rmax.mark(fieldBoundaryMarker, rmax_indx) # Mark the upper radial
                                                       # boundary with rmax_indx
        else:
            fieldBoundaryMarker = None

        # Particle boundary conditions

        # Create a MeshFunction object that contains the mesh, and is
        # defined on mesh facets.  The function has a (size_t) value
        # of 0, 1, or 2, to identify the boundary-condition applied.
        # The function values are set later.

        # A boundary has dimension 1 less than the domain:
        particleBoundaryMarker = df_M.MeshFunction('size_t', mesh, mesh.topology().dim()-1)

        # Initialize all mesh facets with a default value of 0
        # Then overwrite boundary facets that have callback functions.
        particleBoundaryMarker.set_all(0)

        if particleBoundaryDict is not None:
            # Create mesh subdomains to apply boundary-conditions
            Gamma_rmin = XBoundary(rmin)  # Gamma_rmin is the lower radial boundary
            Gamma_rmax = XBoundary(rmax)  # Gamma_rmax is the upper radial boundary
            # At this point, y goes from 0 to 1
            Gamma_thmin = YBoundary(0.0)  # Gamma_thmin is the lower theta boundary
            Gamma_thmax = YBoundary(1.0) # Gamma_thmax is the upper theta boundary

            # Mark these subdomains (boundaries) with non-zero integers
            rmin_indx = particleBoundaryDict['rmin']
            rmax_indx = particleBoundaryDict['rmax']
            thmin_indx = particleBoundaryDict['thmin']
            thmax_indx = particleBoundaryDict['thmax']
            Gamma_rmin.mark(particleBoundaryMarker, rmin_indx) # Mark the lower radial
                                                               # boundary with rmin_indx
            Gamma_rmax.mark(particleBoundaryMarker, rmax_indx) # Mark the upper radial
                                                               # boundary with rmax_indx
            Gamma_thmin.mark(particleBoundaryMarker, thmin_indx) # Mark the lower theta
                                                                 # boundary with thmin_indx
            Gamma_thmax.mark(particleBoundaryMarker, thmax_indx) # Mark the lower theta
                                                                 # boundary with thmax_indx

# Now transform the rectangular mesh to its desired shape

        # Pull out the x, y coordinates of the rectangle
        x = mesh.coordinates()[:,0]
        y = mesh.coordinates()[:,1]

# Make the stretched mesh coordinates

        # Make denser mesh towards r=rmin (x is the r-coordinate)
        def denser(x, y):
            return [rmin + (rmax-rmin)*((x-rmin)/(rmax-rmin))**stretch, y]

        x_bar, y_bar = denser(x, y)

        # Transpose to get them 'interleaved' as x1, y1, x2, y2, ... again
        xy_bar_coor = np_M.array([x_bar, y_bar]).transpose()

        # Put the coords back into the mesh to get a plot
        mesh.coordinates()[:] = xy_bar_coor

# Plot the stretched mesh
#        if plot_flag:
#            df_M.plot(mesh, title='stretched mesh', axes=True)
#            df_M.interactive()
#raw_input('B: Press <ENTER> to continue')

        # Convert x, y coordinates to r, theta coordinates
        def cylinder(r, t):
            return [r*np_M.cos(tmax*t), r*np_M.sin(tmax*t)]

        x_hat, y_hat = cylinder(x_bar, y_bar)

        # Put the coords back into the mesh object to get a plot
        xy_hat_coor = np_M.array([x_hat, y_hat]).transpose()
        mesh.coordinates()[:] = xy_hat_coor

        # Make a plot of the mesh, with non-zero values showing marked
        # boundaries
        if (plot_flag):
            df_M.plot(mesh, title='cylindrical mesh', axes=True)
            df_M.plot(fieldBoundaryMarker, title='field boundary marks', axes=True)
            df_M.plot(particleBoundaryMarker, title='particle boundary marks', axes=True)
            df_M.interactive()
#raw_input('C: Press <ENTER> to continue')

        # Save the class attributes
        self.mesh = mesh
        self.field_boundary_marker = fieldBoundaryMarker
        self.particle_boundary_marker = particleBoundaryMarker

#A copy of a mesh may be created as follows:
#mesh_copy = Mesh(mesh)

        return
#    def create_mesh(self, meshInputCI, plot_flag):ENDDEF

#class UserMesh_C(Mesh_C): ENDCLASS

# The field-solve class could be in a different module file.  Here,
# the connection between the field mesh and the field solve is close, so
# they're both in the same file.

#
# Solve the equations for the fields
#
class UserPoissonSolve_C(PoissonSolve_C):
    """UserPoissonSolve_C solves Poisson's equation to obtain the electric field from
       the charge density.  It can be edited by the user to specify the field
       solvers and any tuning parameters.  The units are MKS by
       default (i.e., if no conversion factor is applied), but any
       units can be used provided the conversion to MKS is available
       in the UserUnits_M.py module.
    """

# Select the unit system to be used for input parameters.
    Convert = U_M.MyPlasmaUnits_C

    def __init__(self, phi, linear_solver, preconditioner, fieldBoundaryMarker, phi_BCs, chargeDensity=None, negElectricField=None):
#    def __init__(self, phi, linear_solver, preconditioner, boundary_marker, phi_rmin, phi_rmax, chargeDensity=None, negElectricField=None):
#    def __init__(self, function_space, linear_solver, preconditioner, boundary_marker, phi_rmin, phi_rmax, computeEflag, mesh=None):
        """Constructor for a Poisson solver written by the user.
        """

        print "This is DOLFIN Version", df_M.DOLFIN_VERSION_STRING

        self.u = phi.function
        V = phi.function_space

#        u_rmin = df_M.Constant(phi_rmin)
#        u_rmax = df_M.Constant(phi_rmax)

        self.charge_density = chargeDensity
        self.neg_electric_field = negElectricField

        # Field-solver parameters
        self.solver_parameters = {}
        if linear_solver is not None:
            self.solver_parameters['linear_solver'] = linear_solver
        if preconditioner is not None:
            self.solver_parameters['preconditioner'] = preconditioner

# Create the Dirichlet boundary-condition list for the Laplace PDE
#       args: DirichletBC(FunctionSpace, GenericFunction, MeshFunction, int, method="topological")
#        V = function_space

        # Get the Dirichlet boundary indices and values
        (rminIndx, phi_rmin) = phi_BCs['rmin']
        (rmaxIndx, phi_rmax) = phi_BCs['rmax']
        
        # Create a function from the boundary values
        uRmin = df_M.Constant(phi_rmin)
        uRmax = df_M.Constant(phi_rmax)

        self.bcs = [df_M.DirichletBC(V, uRmin, fieldBoundaryMarker, rminIndx), df_M.DirichletBC(V, uRmax, fieldBoundaryMarker, rmaxIndx)]

# Define the variational problem
        w = df_M.TrialFunction(V)
        self.v = df_M.TestFunction(V)
        f = df_M.Constant(0.0)
        self.a = df_M.inner(df_M.grad(w), df_M.grad(self.v))*df_M.dx

        df_M.set_log_level(df_M.PROGRESS) # Gives PETSc LU solver, (null). (Direct solver).
#df.set_log_level(1) # Gives the most messages

# default LU is flakey: different answers on different calls: NO!: this was a heap problem of unitialized memory! Fix was to used set_all(0) above.
#        self.phi = None
#        self.negE = None
        return
#    def __init__(self, phi, linear_solver, preconditioner, field_boundary_marker, phi_rmin, phi_rmax, chargeDensity=None, negElectricField=None):ENDDEF

#class UserPoissonSolve_C(PoissonSolve_C):ENDCLASS
