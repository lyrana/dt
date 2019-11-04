# Make a mesh that DT can use

#__version__ = 0.1
#__author__ = 'Copyright (C) 2016 L. D. Hughes'
__all__ = ['UserMesh2DCirc_C',
           'UserMeshInput2DCirc_C', 
           'UserPoissonSolve2DCirc_C',           
#           'UserPoissonSolveInput2DCirc_C',
           ]

"""UserMesh_y_Fields defines the mesh, marks boundaries for field and
   particle boundary-conditions (BCs), and creates the field solvers.
   The type of mesh may be FE or FD.  In this case, it's easiest to
   apply the BCs while the mesh is being created, since at first, the
   boundaries are Cartesian coordinate lines that are easy to specify.
   The mesh is then stretched and transformed to a circular shape.
"""

import sys
import os
import numpy as np_m
import math
# !!! Direct invocation of dolfin. OK because UserMesh2DCirc is a
# sub-class of Mesh_C !!!
import dolfin as df_m
import matplotlib.pyplot as mplot_m

from Dolfin_Module import Mesh_C
from Dolfin_Module import PoissonSolve_C

import UserUnits_Module as U_M

#STARTCLASS
class UserMeshInput2DCirc_C(object):
#class DTmeshInput_C(object):
    """Input for the field mesh
       The user can modify this for different mesh specifications.
       Field mesh, field solve control?  Could use to pass control things to the field mesh
    """

    def __init__(self):
        """ List the mesh variables that the user can set in MAIN.py
        """
#        self.mesh_type_options = ['FE', 'Cartesian']
#        self.mesh_type = None

        self.mesh_file = None

        self.user_mesh_input = None
        self.user_mesh_class = None

        self.precision = None
        self.mesh_class = None

        self.rmin = None
        self.rmax = None
        self.nr = None
        self.stretch = None
        
        self.tmax = None
        self.nt = None

        # Options: 'left, 'right', 'left/right', 'crossed'
        self.diagonal = None

        self.field_boundary_file = None
        # User-assigned names of mesh boundaries where Dirichlet
        # values are set.
        self.field_boundary_dict = None

        self.particle_boundary_file = None
        # User-assigned names of mesh boundaries where particle BCs
        # are set.
        self.particle_boundary_dict = None

        self.particle_source_file = None
        # User-assigned names of mesh regions where particles are
        # created
        self.particle_source_dict = None

# May want things like this in order to call DT from a loop?
# or spawn off many runs?
# maybe don't need all of these:
        self.mesh_M = None
        self.pmesh_M = None

        # the particle mesh is a copy of the field mesh
#        self.pmesh_M = df_m.Mesh(mesh_M)

        return

# Create definitions of functions to test if points are on the boundaries
# The XBoundary function takes care of the two Dirichlet boundaries, and
# set_all() marks every facet with a 2.

# Don't actually need the AllFacets function for this problem
class AllFacets(df_m.SubDomain):
    """The AllFacets class is a specialized SubDomain
    """
    # Return true for every location x that is on a boundary
    def inside(self, x, on_boundary):
        return on_boundary

class XBoundary(df_m.SubDomain):
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

class YBoundary(df_m.SubDomain):
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

#class UserMeshInput2DCirc_C(object):ENDCLASS

# User exposes whatever mesh parameters are useful in __init__ and
# these can be set in __main__

#STARTCLASS
class UserMesh2DCirc_C(Mesh_C):
    """UserMesh2DCirc_C is derived from Mesh_C, and includes marked
       subdomains (e.g., Dirichlet boundaries).  It is to be edited by
       the user to specify the simulation mesh.

       The units are MKS by default (i.e., if no conversion
       factor is applied), but any units can be used provided the
       conversion to MKS is available in the UserUnits_M.py module.
    """

    # Select the unit system to be used for input parameters.
#    Convert = U_M.MyPlasmaUnits_C

    def __init__(self, mesh_input, compute_dictionaries=False, compute_tree=False, plot_flag=False, plot_title=None):
        """
            The class UserMesh2DCirc contains these attributes:
                1. A mesh.
                2. A fieldBoundaryMarker (MeshFunction) that marks
                   mesh boundaries for field boundary-conditions.
                3. A particleBoundaryMarker (MeshFunction) that marks
                   mesh boundaries for particle boundary-conditions.
        """

        coordinateSystem = 'Cartesian' # Mesh coordinates are Cartesian x-y, even though
                                       # the mesh is a 2D quarter-circle.

        if mesh_input.mesh_file is None:
            self.create_mesh(mesh_input, plot_flag=plot_flag, plot_title=plot_title)
            # don't need another mesh plot
            plot_flag = False

        # Call the parent constructor to complete setting class variables.
        mesh_file = mesh_input.mesh_file
#        super(self.__class__, self).__init__(mesh_file=mesh_file, coordinate_system=coordinateSystem, compute_dictionaries=compute_dictionaries, compute_tree=compute_tree, plot_flag=plot_flag)
        super(self.__class__, self).__init__(mesh_file=mesh_file, coordinate_system=coordinateSystem, compute_dictionaries=compute_dictionaries, compute_tree=compute_tree)
        
        self.field_boundary_dict = mesh_input.field_boundary_dict
        self.particle_boundary_dict = mesh_input.particle_boundary_dict

        # Read boundary markers from file, if provided
        fieldBoundaryFile = mesh_input.field_boundary_file
        particleBoundaryFile = mesh_input.particle_boundary_file
        if mesh_input.field_boundary_file is not None:
            fieldBoundaryMarker = df_m.MeshFunction("size_t", self.mesh, fieldBoundaryFile)
            self.field_boundary_marker = fieldBoundaryMarker
        if particleBoundaryFile is not None:
            particleBoundaryMarker = df_m.MeshFunction("size_t", self.mesh, particleBoundaryFile)
# or equivalently:       particleBoundaryMarker = df_m.MeshFunctionSizet(self.mesh, "Pbcs_quarter_circle_mesh_crossed.xml")
            self.particle_boundary_marker = particleBoundaryMarker

        return  
#    def __init__(self, mesh_input=None, compute_dictionaries=False, compute_tree=False, plot_flag=False):ENDDEF

# Inherited from Mesh_C:
#    def copy(self):
#        return copy.deepcopy(self)

#class UserMesh2DCirc(Mesh_C):
    def create_mesh(self, mesh_input, plot_flag=False, plot_title=None):
        """Create a mesh and mark subdomains (e.g., Dirichlet
           boundaries and particle boundaries) according to the user's
           specifications.
        """

        # Radial limits of mesh
        rmin = mesh_input.rmin
        rmax = mesh_input.rmax

        stretch = mesh_input.stretch
        nr = mesh_input.nr

        # The final mesh will go from 0 to tmax:
        tmax = mesh_input.tmax
        nt = mesh_input.nt

        # Options: 'left, 'right', 'left/right', 'crossed'
        diagonal = mesh_input.diagonal

        # Boundary conditions for fields and particles
        fieldBoundaryDict = mesh_input.field_boundary_dict
        particleBoundaryDict = mesh_input.particle_boundary_dict

        plotTitle = plot_title

# First, make a rectangular mesh

#       args: RectangleMesh(x0, y0, x1, y1, nx, ny, diagonal="right")

# Switch these for CEE/Laptop
#        if df_m.DOLFIN_VERSION_STRING > "1.5.0":
        if df_m.__version__ > "1.5.0":
# v > 1.5:
            mesh = df_m.RectangleMesh(df_m.Point(rmin, 0.0), df_m.Point(rmax, 1.0), nr, nt, diagonal)
        else:
# v = 1.5:
            mesh = df_m.RectangleMesh(rmin, 0.0, rmax, 1.0, nr, nt, diagonal)

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
#            fieldBoundaryMarker = df_m.FacetFunction('size_t', mesh)
            fieldBoundaryMarker = df_m.MeshFunction('size_t', mesh, mesh.topology().dim()-1)

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
        particleBoundaryMarker = df_m.MeshFunction('size_t', mesh, mesh.topology().dim()-1)

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
        xy_bar_coor = np_m.array([x_bar, y_bar]).transpose()

        # Put the coords back into the mesh to get a plot
        mesh.coordinates()[:] = xy_bar_coor

# Plot the stretched mesh
#        if plot_flag:
#            df_m.plot(mesh, title='stretched mesh', axes=True)
#            df_m.interactive()
#raw_input('B: Press <ENTER> to continue')

        # Convert x, y coordinates to r, theta coordinates
        def cylinder(r, t):
            return [r*np_m.cos(tmax*t), r*np_m.sin(tmax*t)]

        x_hat, y_hat = cylinder(x_bar, y_bar)

        # Put the coords back into the mesh object to get a plot
        xy_hat_coor = np_m.array([x_hat, y_hat]).transpose()
        mesh.coordinates()[:] = xy_hat_coor

        # Make a plot of the mesh, with non-zero values showing marked
        # boundaries
        if (plot_flag):
            fileName = os.path.basename(__file__)
            if plot_title is None: plotTitle = fileName + ": RZ-mesh"
#            df_m.plot(mesh, title=plotTitle, axes=True)
            df_m.plot(mesh, title=plotTitle)
#*** Warning: Matplotlib plotting backend does not support mesh function of dim 1. Continuing without plotting...            
#            df_m.plot(fieldBoundaryMarker, title=fileName + ': field boundary marks', axes=True)
#            df_m.plot(particleBoundaryMarker, title=fileName + ': particle boundary marks', axes=True)
#            df_m.interactive()
            mplot_m.show()            
#raw_input('C: Press <ENTER> to continue')

        # Save the class attributes
        self.mesh = mesh
        self.field_boundary_marker = fieldBoundaryMarker
        self.particle_boundary_marker = particleBoundaryMarker

#A copy of a mesh may be created as follows:
#mesh_copy = Mesh(mesh)

        return
#    def create_mesh(self, mesh_input, plot_flag):ENDDEF

#class UserMesh2DCirc(Mesh_C): ENDCLASS

# The field-solve class could be in a different module file.  Here,
# the connection between the field mesh and the field solve is close, so
# they're both in the same file.

#
# Solve the equations for the fields
#
#STARTCLASS
class UserPoissonSolve2DCirc_C(PoissonSolve_C):
    """UserPoissonSolve2DCirc_C solves Poisson's equation on a 2D region that is a
       segment of a circle to obtain the electric field from the charge density.  It
       can be edited by the user to specify the field solvers and any tuning
       parameters.  The units are MKS by default (i.e., if no conversion factor is
       applied), but any units can be used provided the conversion to MKS is
       available in the UserUnits_M.py module.
    """

    def __init__(self, phi_F, linear_solver, preconditioner, field_boundary_marker, phi_BCs, neg_electric_field=None):
        """Constructor for a Poisson solver written by the user.

           :cvar bool pde_has_constant_coeffs: If true, the FE matrix has to assembled
                                               only at initialization.

        """

        # Select the unit system to be used for input parameters.
        constants = U_M.MyPlasmaUnits_C

        # Call the PoissonSolve_C base class constructor for
        # non-problem-specific initialization.
        super(self.__class__, self).__init__(phi_F)

        # Field-solver parameters
        self.solver_parameters = {}
        if linear_solver is not None:
            self.solver_parameters['linear_solver'] = linear_solver
        if preconditioner is not None:
            self.solver_parameters['preconditioner'] = preconditioner

        self.neg_electric_field = neg_electric_field
        
        # Get the Dirichlet boundary indices and values
        (rminIndx, phi_rmin) = phi_BCs['rmin']
        (rmaxIndx, phi_rmax) = phi_BCs['rmax']

        # Create a function from the boundary values
        uRmin = df_m.Constant(phi_rmin)
        uRmax = df_m.Constant(phi_rmax)

        # Create the Dirichlet boundary-condition list for the Laplace PDE
        # args: DirichletBC(FunctionSpace, GenericFunction, MeshFunction, int, method="topological")
        
        self.bcs = [df_m.DirichletBC(self.V, uRmin, field_boundary_marker, rminIndx), df_m.DirichletBC(self.V, uRmax, field_boundary_marker, rmaxIndx)]

        # Define the variational problem

        ## Make the bilinear form 'a(w,v)' ##

        # The bilinear form 'a' is the LHS in the variational form of the
        # Laplace/Poisson eq. It is the 'inner product' of grad(w) times
        # grad(v), i.e., the integral of this product (times any coefficients)
        # over the domain. The spherical-coordinate form has a coefficient r**2
        # (for the volume-element needed to integrate over the domain). This
        # Cartesian-coordinate form doesn't.
        epsilon = df_m.Constant(constants.epsilon_0)
        self.a = epsilon*df_m.inner(df_m.grad(self.w), df_m.grad(self.v))*df_m.dx

        # Specify whether 'a(u,v)' has time-independent coefficients.  If so, the
        # matrix A only needs to be assembled once, avoiding redundant work.
        # By inspection of the expression above:
        self.pde_has_constant_coeffs = True

        # Carry out the initial assembly
        self.assemble_matrix()

        # Initialize the source vector 'b' to zero charge-density.
        self.assemble_source_expression(0.0)

        # Set message detail level for the field-solver        
        df_m.set_log_level(df_m.LogLevel.PROGRESS) # Gives PETSc LU solver, (null). (Direct solver).
#        df_m.set_log_level(16) # Gives PETSc LU solver, (null). (Direct solver).
#df.set_log_level(1) # Gives the most messages

# default LU is flakey: different answers on different calls: NO!: this was a heap problem of unitialized memory! Fix was to used set_all(0) above.
#        self.phi = None
#        self.negE = None

        return
#    def __init__(self, phi, linear_solver, preconditioner, field_boundary_marker, phi_BCs, charge_density=None, neg_electric_field=None):ENDDEF

#class UserPoissonSolve2DCirc_C(PoissonSolve_C):ENDCLASS
