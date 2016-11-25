# Make a mesh that DT can use

#__version__ = 0.1
#__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

"""UserFields defines the mesh, field boundary conditions (BCs), and
   the field solver.  The type of mesh may be FE or FD.  In this case,
   it's easiest to apply the BCs while the mesh is being created,
   since at first, the boundaries are Cartesian coordinate lines that
   are easy to specify.
"""

import sys
import os
import dolfin as df_M
import numpy as np_M
import math

from Dolfin_Module import Mesh_C
from Dolfin_Module import PoissonSolve_C

import UserUnits_Module as U_M

# Define subdomain classes to test if points are on the boundaries
# The XBoundary function takes care of the two Dirichlet boundaries, and
# set_all() marks every facet with a 2.

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
        # on_boundary is passed in with value True for all outer boundaries
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

       This mesh is 1D.

    """

# Select the unit system to be used for input parameters.
    Convert = U_M.MyPlasmaUnits_C

    def __init__(self, meshInputCI=None, mesh=None, boundaryMarker=None, computeDictionaries=False, computeTree=False, plotFlag=False):
        """
            The class UserMesh_C contains these attributes:
                1. A mesh.
                2. A boundaryMarker (MeshFunction) that marks mesh boundaries for boundary conditions.
        """

        if meshInputCI is not None:
            self.create_mesh(meshInputCI, plotFlag)
            # Don't need another mesh plot
            plotFlag = False
        else:
            error_msg = "UserMesh_C: no mesh specified"
            sys.exit(error_msg)

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

# First, make a 1-D mesh
        mesh = df_M.UnitIntervalMesh(nr)

#
# Mark the Dirichlet-boundaries before transforming the mesh.
#

# Create a MeshFunction object that is defined on mesh facets (which
# in this 1D case are cell vertices.)  The function has (size_t) value
# of 0, 1, or 2 here.  These mark the mesh vertices where Dirichlet
# conditions need to be set.  The boundary potentials corresponding to
# these numbers are set later.

        # A 'boundary' by definition has dimension 1 less than the domain:
        boundaryMarker = df_M.MeshFunction("size_t", mesh, mesh.topology().dim()-1)

#        print 'dim =', mesh.topology().dim()

        # Initialize all the facets with a default value of 2.
        # Then overwrite boundary facets that have Dirichlet BCs.
        boundaryMarker.set_all(2)

        # Create mesh subdomains Gamma_nnn that are the boundaries
        # (Could also do this later when the mesh has been stretched)
        Gamma_rmin = XBoundary(0) # Gamma_rmin is the lower radial boundary
        Gamma_rmax = XBoundary(1) # Gamma_rmax is the upper radial boundary

        # Apply the  boundaryMarker function to these subdomains.
        Gamma_rmin.mark(boundaryMarker, 0) # Mark the lower radial boundary with a 0
        Gamma_rmax.mark(boundaryMarker, 1) # Mark the upper radial boundary with a 1
#
# Now transform the unit interval mesh to its desired shape
#
        # Get the x coordinates of the unit interval
        x = mesh.coordinates()[:,0]

#        print 'x =', x
        # Make the stretched mesh coordinates

        # Function to stretch the mesh as r increases from rmin to rmax
        def map_coordinate(x):
            """Maps the [0, 1] interval to [rmin, rmax] with a stretch factor.
            """
            return rmin + (rmax-rmin)*x**stretch

        r = map_coordinate(x)

#        print 'r = ', r

        # Put the coords back into the mesh.
        mesh.coordinates()[:,0] = r

# Make a plot of the mesh
        if (plotFlag):
            df_M.plot(mesh, title='radial mesh', axes=True)
            df_M.interactive()

        self.boundary_marker = boundaryMarker
        self.mesh = mesh

        return


# The field-solve class could be in a different module file.  Here,
# the connection between the field mesh and the field solve is close, so
# they're both in the same file.

#
# Solve the equations for the fields
#
#ClassClassClassClassClassClassClassClassclass
class UserPoissonSolve_C(PoissonSolve_C):
    """UserPoissonSolve_C solves equations to obtain physical fields from
       sources.  It can be edited by the user to specify the field
       solvers and any tuning parameters.  The units are MKS by
       default (i.e., if no conversion factor is applied), but any
       units can be used provided the conversion to MKS is available
       in the UserUnits_M.py module.
    """

# Select the unit system to be used for input parameters.
    Convert = U_M.MyPlasmaUnits_C

    def __init__(self, phi, linear_solver, preconditioner, boundary_marker, phi_rmin, phi_rmax, chargeDensity=None, negElectricField=None):
        """mesh argument is only needed if, e.g., using a SpatialCoordinate in the equations.
        """

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
        self.bcs = [df_M.DirichletBC(V, u_rmin, boundary_marker, 0), df_M.DirichletBC(V, u_rmax, boundary_marker, 1)]

        # Define the variational problem
        w = df_M.TrialFunction(V)
        self.v = df_M.TestFunction(V)
        f = df_M.Constant(0.0)
        r = df_M.SpatialCoordinate(V.mesh())

        # Spherical-coordinate form has factor r**2
        self.a = df_M.inner(df_M.nabla_grad(w), df_M.nabla_grad(self.v))*r[0]**2*df_M.dx

        df_M.set_log_level(df_M.PROGRESS) # Gives PETSc LU solver, (null). (Direct solver).
#df.set_log_level(1) # Gives the most messages

# default LU is flakey: different answers on different calls: NO!: this was a heap problem of unitialized memory!
#        self.phi = None
#        self.negE = None

        return
