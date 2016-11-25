# Make a mesh that DT can use

#__version__ = 0.1
#__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

"""UserFields defines the mesh, field boundary conditions (BCs), and
   the field solvers.  The type of mesh may be FE or FD.  In this case,
   it's easiest to apply the BCs while the mesh is being created,
   since at first, the boundaries are Cartesian coordinate lines that
   are easy to specify.  The mesh is then stretched and transformed to
   a cylindrical shape.
"""

import sys
import os
import dolfin as df_M
import numpy as np_M
import math

from Dolfin_Module import Mesh_C
from Dolfin_Module import PoissonSolve_C

import UserUnits_Module as U_M

# Create definitions of functions to test if points are on the boundaries
# The XBoundary function takes care of the two Dirichlet boundaries, and
# set_all() marks every facet with a 2.

# Don't actually need the AllFacets function for this problem
#ClassClassClassClassClassClassClassClassclass
class AllFacets(df_M.SubDomain):
    """The AllFacets class is a specialized SubDomain
    """
    # Return true for every location x that is on a boundary
    def inside(self, x, on_boundary):
        return on_boundary

#ClassClassClassClassClassClassClassClassclass
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


# User exposes whatever mesh parameters are useful in __init__ and
# these can be set in __main__

#ClassClassClassClassClassClassClassClassclass
class UserMesh_C(Mesh_C):
    """UserMesh_C is derived from Mesh_C.  It is to be edited by the user to specify the simulation
       mesh.  The units are MKS by default (i.e., if no conversion
       factor is applied), but any units can be used provided the
       conversion to MKS is available in the UserUnits_M.py module.
    """

    # Select the unit system to be used for input parameters.
    Convert = U_M.MyPlasmaUnits_C

    def __init__(self, meshInputCI=None, mesh=None, boundaryMarker=None, computeDictionaries=False, computeTree=False, plotFlag=False):
        """
            The class UserMesh_C contains these attributes:
                1. A mesh.
                2. A MeshFunction that marks mesh boundaries for boundary conditions.
        """

        if meshInputCI is not None:
            self.create_mesh(meshInputCI, plotFlag)
            # don't need another mesh plot
            plotFlag = False
        elif mesh is not None:
            self.mesh = mesh
            self.boundary_marker = boundaryMarker
        else:
            error_msg = "UserMesh_C: no mesh specified"
            sys.exit(error_msg)

#        self.gdim = self.mesh.geometry().dim()
#        self.tdim = self.mesh.topology().dim()
        # Call the parent constructor to complete setting class variables.
        super(self.__class__, self).__init__(meshFile=None, computeDictionaries=computeDictionaries, computeTree=computeTree, plotFlag=plotFlag)

        return

    

    def create_mesh(self, miCI, plotFlag):
        """
           Create a mesh according to the user's specifications.
        """
        rmin = miCI.rmin
        rmax = miCI.rmax
        stretch = miCI.stretch
        nr = miCI.nr

        tmax = miCI.tmax
        nt = miCI.nt

        # Options: 'left, 'right', 'left/right', 'crossed'
        diagonal = miCI.diagonal
        
# The final mesh will go from 0 to tmax:
#    Theta = math.pi/2 # theta extent of the mesh
        Theta = tmax # theta extent of the mesh

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

#
# Set the Dirichlet boundary conditions on the Rectangle, before transforming it,
# since it's a simple rectangle then.
#

# Create a MeshFunction object that contains the mesh, and is defined on mesh facets.
# The function has (size_t) value of 0, 1, or 2, to identify the boundary-condition applied.
# The function values are set later.

        # A boundary has dimension 1 less than the domain:
        boundaryMarker = df_M.MeshFunction("size_t", mesh, mesh.topology().dim()-1)

        # Initialize all mesh facets with a default value of 2
        # Then overwrite boundary facets that have Dirichlet BCs.
        boundaryMarker.set_all(2)

# Create the functions used to mark boundary facets by instantiating
# an XBoundary object
        Gamma_rmin = XBoundary(rmin) # Gamma_rmin is the lower radial boundary
        Gamma_rmax = XBoundary(rmax) # Gamma_rmax is the upper radial boundary

        # Use the functions to mark the mesh
        Gamma_rmin.mark(boundaryMarker, 0) # Mark the lower radial boundary with a 0
        Gamma_rmax.mark(boundaryMarker, 1) # Mark the upper radial boundary with a 1

#
# Now transform the rectangular mesh to its desired shape
#

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
#        if plotFlag:
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

        # Make a plot of the mesh
        if (plotFlag):
            df_M.plot(mesh, title='cylindrical mesh', axes=True)
            df_M.interactive()
#raw_input('C: Press <ENTER> to continue')

        self.boundary_marker = boundaryMarker
        self.mesh = mesh
#        self.mesh = mesh_type

#A copy of a mesh may be created as follows:

#mesh_copy = Mesh(mesh)

# Compute the search tree. Needed if there are field probes,
# etc. where fields need to be evaluated at arbitrary positions within
# the mesh.

# The particle mesh needs a tree to compute the forces on the
# particles.

        # compute_tree = True

        # if (compute_tree):
        #     self.mesh_tree = self.mesh.bounding_box_tree()

            # endif


        return


# The field-solve class could be in a different module file.  Here,
# the connection between the field mesh and the field solve is close, so
# they're both in the same file.

#
# Solve the equations for the fields
#
#ClassClassClassClassClassClassClassClassclass
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

    def __init__(self, phi, linear_solver, preconditioner, boundary_marker, phi_rmin, phi_rmax, chargeDensity=None, negElectricField=None):
#    def __init__(self, function_space, linear_solver, preconditioner, boundary_marker, phi_rmin, phi_rmax, computeEflag, mesh=None):
        """Constructor for a Poisson solver written by the user.
        """

        print "This is DOLFIN Version", df_M.DOLFIN_VERSION_STRING

        self.u = phi.function
        V = phi.function_space

        u_rmin = df_M.Constant(phi_rmin)
        u_rmax = df_M.Constant(phi_rmax)

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
        self.bcs = [df_M.DirichletBC(V, u_rmin, boundary_marker, 0), df_M.DirichletBC(V, u_rmax, boundary_marker, 1)]

# Define the variational problem
        w = df_M.TrialFunction(V)
        self.v = df_M.TestFunction(V)
        f = df_M.Constant(0.0)
        self.a = df_M.inner(df_M.grad(w), df_M.grad(self.v))*df_M.dx

        df_M.set_log_level(df_M.PROGRESS) # Gives PETSc LU solver, (null). (Direct solver).
#df.set_log_level(1) # Gives the most messages

# default LU is flakey: different answers on different calls: NO!: this was a heap problem of unitialized memory!
#        self.phi = None
#        self.negE = None

        return
