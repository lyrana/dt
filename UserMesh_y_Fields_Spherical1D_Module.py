# Make a mesh that DT can use

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
__all__ = ['UserMesh1DS_C',
           'UserMeshInput1DS_C', 
           'UserPoissonSolve1DS_C',           
           'UserPoissonSolveInput1DS_C',
           ]



"""UserFields defines the mesh, field boundary conditions (BCs), and
   the field solver.  The type of mesh may be FE or FD.  In this case,
   it's easiest to apply the BCs while the mesh is being created,
   since at first, the boundaries are Cartesian coordinate lines that
   are easy to specify.
"""

import sys
import os
#import numpy as np_m
import math
# !!! Direct invocation of dolfin. OK because UserMesh1DS_C is a
# sub-class of Mesh_C !!!
import dolfin as df_m
import matplotlib.pyplot as mplot_m

from Dolfin_Module import Mesh_C
from Dolfin_Module import PoissonSolve_C

import UserUnits_Module as U_M

#STARTCLASS
class UserMeshInput1DS_C(object):
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
        
        self.field_boundary_file = None
        # User-assigned names of mesh boundaries where Dirichlet
        # values are set.
        self.field_boundary_dict = None
        # The mesh function that marks facets with BCs
        self.field_boundary_marker = None

        self.particle_boundary_file = None
        # User-assigned names of mesh boundaries where particle BCs
        # are set.
        self.particle_boundary_dict = None
        # The mesh function that marks facets with BCs
        self.particle_boundary_marker = None

        self.particle_source_file = None
        # User-assigned names of mesh regions where particles are
        # created
        self.particle_source_dict = None
        # The mesh function that marks source cells
        self.particle_source_marker = None

# May want things like this in order to call DT from a loop?
# or spawn off many runs?
# maybe don't need all of these:
        self.mesh_M = None
        self.pmesh_M = None

        # the particle mesh is a copy of the field mesh
#        self.pmesh_M = df_m.Mesh(mesh_M)

        return
        
#class UserMeshInput1DS_C(object):
    def __str__(self):
        """Print the class members.
        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():'
        print fncName

        print_string = ' rmin = ' + str(self.rmin)
        print_string += '\n rmax = ' + str(self.rmax)
        print_string += '\n stretch = ' + str(self.stretch)
        print_string += '\n nr = ' + str(self.nr)
        print_string += '\n field_boundary_dict = ' + str(self.field_boundary_dict)
        print_string += '\n particle_boundary_dict = ' + str(self.particle_boundary_dict)
        print_string += '\n particle_boundary_file = ' + str(self.particle_boundary_file)

        return print_string
#     def __str__(self):ENDDEF

#class UserMeshInput1DS_C(object):ENDCLASS


# Define subdomain classes to test if points are on the boundaries
# The XBoundary function takes care of the two Dirichlet boundaries, and
# set_all() marks every facet with a 2.

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
        # on_boundary is passed in with value True for all outer boundaries
        tol = 1.0e-10
        return on_boundary and abs(x[0]-self.x_boundary_value) < tol

#class XBoundary(df_m.SubDomain):ENDCLASS

# User exposes whatever mesh parameters are useful in __init__ and
# these can be set in __main__

#STARTCLASS
class UserMesh1DS_C(Mesh_C):
    """UserMesh1DS_C is derived from Mesh_C.  It is to be edited by the user to specify the simulation
       mesh.  The units are MKS by default (i.e., if no conversion
       factor is applied), but any units can be used provided the
       conversion to MKS is available in the UserUnits_M.py module.

       This mesh is 1D.

    """

# Select the unit system to be used for input parameters.
    Convert = U_M.MyPlasmaUnits_C

# Mesh_C constructor:
#class UserMesh1DS_C(Mesh_C):
    def __init__(self, mesh_input, compute_dictionaries=False, compute_tree=False, plot_flag=False, plot_title=None):
        """
            The class UserMesh1DS_C contains these attributes:
                1. A mesh.
                2. A fieldBoundaryMarker (MeshFunction) that marks
                   mesh boundaries for field boundary-conditions.
                3. A particleBoundaryMarker (MeshFunction) that marks
                   mesh boundaries for particle boundary-conditions.
        """

        if mesh_input.mesh_file is None:
            self.create_mesh(mesh_input, plot_flag=plot_flag, plot_title=plot_title)
            # Don't need another mesh plot
            plot_flag = False

        # Call the parent constructor to complete setting class variables.
        meshFile = mesh_input.mesh_file
        super(self.__class__, self).__init__(mesh_file=meshFile, compute_dictionaries=compute_dictionaries, compute_tree=compute_tree, plot_flag=plot_flag, plot_title=plot_title)

        self.field_boundary_dict = mesh_input.field_boundary_dict
        self.particle_boundary_dict = mesh_input.particle_boundary_dict

        # Read boundary markers from file, if provided
        fieldBoundaryFile = mesh_input.field_boundary_file
        particleBoundaryFile = mesh_input.particle_boundary_file
        if fieldBoundaryFile is not None:
            fieldBoundaryMarker = df_m.MeshFunctionSizet(self.mesh, fieldBoundaryFile)
        if particleBoundaryFile is not None:
            particleBoundaryMarker = df_m.MeshFunctionSizet(self.mesh, particleBoundaryFile)
# or:       particleBoundaryMarker = df_m.MeshFunctionSizet(self.mesh, "Pbcs_quarter_circle_mesh_crossed.xml")

        return
#    def __init__(self, mesh_input, compute_dictionaries=False, compute_tree=False, plot_flag=False):ENDDEF

#class UserMesh1DS_C(Mesh_C):
    def __str__(self):
        """Print the class members.
        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():'
        print fncName

        print_string = ' mesh_tdim = ' + str(self.mesh_tdim)
        print_string += '\nmesh_gdim = ' + str(self.mesh_gdim)

        return print_string
#     def __str__(self):ENDDEF


# Inherited from Mesh_C:
#    def copy(self):
#        return copy.deepcopy(self)

#class UserMesh1DS_C(Mesh_C):
    def create_mesh(self, mesh_input, plot_flag=False, plot_title=None):
        """
           Create a mesh according to the user's specifications.
        """

        rmin = mesh_input.rmin
        rmax = mesh_input.rmax

        fieldBoundaryDict = mesh_input.field_boundary_dict
        particleBoundaryDict = mesh_input.particle_boundary_dict

        stretch = mesh_input.stretch
        nr = mesh_input.nr

        plotTitle = plot_title

        # First, make a 1-D mesh
        mesh = df_m.UnitIntervalMesh(nr)

        # Field boundary conditions

        if fieldBoundaryDict is not None:
            #
            # Mark the Dirichlet-boundaries before transforming the mesh.
            #

            # Create a MeshFunction object that is defined on mesh facets (which
            # in this 1D case are cell vertices.)  The function has (size_t) value
            # of 0, 1, or 2 here.  These mark the mesh vertices where Dirichlet
            # conditions need to be set.  The boundary potentials corresponding to
            # these numbers are set later.

            # A 'boundary' by definition has dimension 1 less than the domain:
            fieldBoundaryMarker = df_m.MeshFunction('size_t', mesh, mesh.topology().dim()-1)

#        print 'dim =', mesh.topology().dim()

            # Initialize all the facets with a default value of 0.
            # Then overwrite boundary facets that have Dirichlet BCs.
            fieldBoundaryMarker.set_all(0)

            # Create mesh subdomains Gamma_nnn that are the boundaries
            # (Could also do this later when the mesh has been stretched)
            # We used UnitInterval to create the mesh, so the boundaries
            # are at 0.0 and 1.0.
            Gamma_rmin = XBoundary(0.0) # Gamma_rmin is the lower radial boundary
            Gamma_rmax = XBoundary(1.0) # Gamma_rmax is the upper radial boundary

            # Assign integer IDs to the boundaries with Dirichlet
            # conditions.  Create a FacetFunction to store these values.
            rmin_indx = fieldBoundaryDict['rmin']
            rmax_indx = fieldBoundaryDict['rmax']
            Gamma_rmin.mark(fieldBoundaryMarker, rmin_indx) # Mark the lower
                                                            # radial boundary
                                                            # with rmin_indx
            Gamma_rmax.mark(fieldBoundaryMarker, rmax_indx) # Mark the upper
                                                            # radial boundary
                                                            # with rmax_indx
        else:
            fieldBoundaryMarker = None
        #   END:if fieldBoundaryDict is not None:


        # Particle boundary conditions are set similarly

        # Mark the facets for particle boundary conditions
        # Create a MeshFunction object defined on the facets of this
        # mesh.  The function has a default (size_t) value of 0,
        # meaning 'no action needed'

        # A boundary has dimension 1 less than the domain:
        particleBoundaryMarker = df_m.MeshFunction('size_t', mesh, mesh.topology().dim()-1)
        # Initialize all mesh facets with a default value of 0.
        particleBoundaryMarker.set_all(0)

        if particleBoundaryDict is not None:
            particleBoundaryMarker = df_m.MeshFunction('size_t', mesh, mesh.topology().dim()-1)
            particleBoundaryMarker.set_all(0)
            Gamma_rmin = XBoundary(0.0) # Gamma_rmin is the lower radial boundary
            Gamma_rmax = XBoundary(1.0) # Gamma_rmax is the upper radial boundary

            rmin_indx = particleBoundaryDict['rmin']
            rmax_indx = particleBoundaryDict['rmax']
            Gamma_rmin.mark(particleBoundaryMarker, rmin_indx)
            Gamma_rmax.mark(particleBoundaryMarker, rmax_indx)
        else:
            particleBoundaryMarker = None
        #END:if particleBoundaryDict is not None:

# Now transform the unit interval mesh to its desired shape

        # Get the x coordinates of the unit interval
        x = mesh.coordinates()[:,0]

#        print 'x =', x
        # Make the stretched mesh coordinates

        # Function to stretch the mesh as r increases from rmin to rmax
        def map_coordinate(x):
            """Maps the [0.0, 1.0] interval to [rmin, rmax] with a stretch factor.
            """
            return rmin + (rmax-rmin)*x**stretch

        r = map_coordinate(x)

#        print 'r = ', r

        # Put the coords back into the mesh.
        mesh.coordinates()[:,0] = r

# Make a plot of the mesh
        if (plot_flag):
            if plot_title is None: plotTitle = "1D radial"
            df_m.plot(mesh, title=plotTitle + " mesh")
            mplot_m.show()

        # Save the class attributes
        self.mesh = mesh
        self.field_boundary_marker = fieldBoundaryMarker
        self.particle_boundary_marker = particleBoundaryMarker

        return
#    def create_mesh(self, mesh_input, plot_flag=False, plot_title=None):ENDDEF

#class UserMesh1DS_C(Mesh_C):ENDCLASS

# The field-solve class could be in a different module file.  Here,
# the connection between the field mesh and the field solve is close, so
# they're both in the same file.

# This isn't used anywhere: not enough args to bother with an input class?

#STARTCLASS
class UserPoissonSolveInput1DS_C(object):
    """Input for the field solver(s).
       The user can modify this for different field solvers.
    """

    def __init__(self):
        """ List the field-solver parameters that the user
        can set in MAIN.py
        """

        self.user_poissonsolve_input = None
        self.user_poissonsolve_class = None

        self.mesh_M = None

        self.element_type = None
        self.element_degree = None

        self.linear_solver = None
        self.preconditioner = None

        # Dirichlet BC object
        self.phi_BCs = None

        self.computeEflag = None

        return

#class UserPoissonSolveInput1DS_C(object):ENDCLASS

#
# Solve the equations for the fields
#
#STARTCLASS
class UserPoissonSolve1DS_C(PoissonSolve_C):
    """UserPoissonSolve1DS_C solves equations to obtain physical fields from
       sources.  It can be edited by the user to specify the field
       solvers and any tuning parameters.  The units are MKS by
       default (i.e., if no conversion factor is applied), but any
       units can be used provided the conversion to MKS is available
       in the UserUnits_M.py module.
    """

    def __init__(self, phi, linear_solver, preconditioner, field_boundary_marker, phi_BCs, charge_density=None, assembled_charge=None, neg_electric_field=None):
        """Initialize a Poisson solver for 1D spherical coordinates

           Generates the bilinear form a(u,v) for the variational form of Poisson's eq.

        """

        # Select the unit system to be used for input parameters.
        constants = U_M.MyPlasmaUnits_C

        self.u = phi.function
        V = phi.function_space

        # Get the Dirichlet boundary values
        (rmin_indx, phi_rmin) = phi_BCs['rmin']
        (rmax_indx, phi_rmax) = phi_BCs['rmax']

        self.charge_density = charge_density
        self.assembled_charge = assembled_charge
        self.neg_electric_field = neg_electric_field

        # Field-solver parameters
        self.solver_parameters = {}
        if linear_solver is not None:
            self.solver_parameters['linear_solver'] = linear_solver
        if preconditioner is not None:
            self.solver_parameters['preconditioner'] = preconditioner

        # Create the Dirichlet boundary-condition list for the Laplace
        # PDE.  The indices rmin_indx and rmax_indx identify the
        # facets in boundary_marker where boundary values are to be
        # set.
        #       args: DirichletBC(FunctionSpace, GenericFunction, MeshFunction, int, method="topological")

        u_rmax = df_m.Constant(phi_rmax)
        if phi_rmin == 'unset':
            # grad phi = 0 at rmin: natural BC.
            self.bcs = [df_m.DirichletBC(V, u_rmax, field_boundary_marker, rmax_indx), ]
        else:
            # phi is set explicitly at both boundaries
            u_rmin = df_m.Constant(phi_rmin)
            self.bcs = [df_m.DirichletBC(V, u_rmin, field_boundary_marker, rmin_indx), df_m.DirichletBC(V, u_rmax, field_boundary_marker, rmax_indx)]
            
        ### Set up the variational problem ###

        w = df_m.TrialFunction(V)
        self.v = df_m.TestFunction(V)

        ## Make the bilinear form 'a(w,v)' ##

        # The bilinear form 'a' is the LHS in the variational form of the
        # Laplace/Poisson eq. It is the 'inner product' of grad(w) times
        # grad(v), i.e., the integral of this product (times any coefficients)
        # over the domain. The spherical-coordinate form has a coefficient r**2
        # (for the volume-element needed to integrate over the domain).

#        f = df_m.Constant(0.0)
        epsilon = df_m.Constant(constants.epsilon_0)
        r = df_m.SpatialCoordinate(V.mesh())

        self.a = epsilon*df_m.inner(df_m.nabla_grad(w), df_m.nabla_grad(self.v))*r[0]**2*df_m.dx

        # Specify whether 'a' has time-independent coefficients.
        self.pde_has_constant_coeffs = True

        # If so, the A matrix only needs to be assembled once.
# Moved
        # if self.pde_has_constant_coeffs is True:
        #     self.A = df_m.assemble(self.a)
        #     print "A is of type", type(self.A)
        #     for bc in self.bcs:
        #         bc.apply(self.A)

        
        #df_m.set_log_level(1) # Gives the most messages
        #df_m.set_log_level(df_m.PROGRESS) # Gives PETSc LU solver, (null). (Direct solver).
        # Turn off solver messages
        df_m.set_log_active(False)

# default LU is flakey: different answers on different calls: NO!: this was a heap problem of unitialized memory!
#        self.phi = None
#        self.negE = None

        # Call the PoissonSolve_C base class constructor for
        # non-problem-specific initialization.
        super(self.__class__, self).__init__()

        return
#    def __init__(self, phi, linear_solver, preconditioner, field_boundary_marker, phi_rmin, phi_rmax, chargeDensity=None, negElectricField=None):ENDDEF

#class UserPoissonSolve1DS_C(PoissonSolve_C):
    def __str__(self):
        """Print the class members.
        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():'
        print fncName

        print_string = 'solver_parameters = ' + str(self.solver_parameters)
        print_string += '\npde_has_constant_coeffs = ' + str(self.pde_has_constant_coeffs)
        print_string += '\ncharge_density = ' + str(self.charge_density)

        return print_string
#    def __str__(self):ENDDEF

#class UserPoissonSolve1DS_C(PoissonSolve_C):ENDCLASS
