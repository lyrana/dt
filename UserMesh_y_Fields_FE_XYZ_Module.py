# Make a mesh that DT can use

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
__all__ = ['UserMesh_C',
           'UserMeshInput_C',
           'UserPoissonSolve1D_C',
           ]

"""UserMesh defines the mesh.
"""

import sys
import os
#import numpy as np_m
import math
# !!! Direct invocation of dolfin. OK because UserMesh_C is a
# sub-class of Mesh_C !!!
import dolfin as df_m
import matplotlib.pyplot as mplot_m

from Dolfin_Module import Mesh_C
from Dolfin_Module import CellSet_C
from Dolfin_Module import PoissonSolve_C

import UserUnits_Module as U_M

### Define boundary geometry

#STARTCLASS
class UserMeshInput_C(object):
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

        self.precision = None

        self.pmin = None
        self.pmax = None

        self.cells_on_side = []

        # 'diagonal' only applies to 2D rectangle:
        # Options: 'left, 'right', 'left/right', 'crossed'
        self.diagonal = None

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

        return

#class UserMeshInput_C(object):ENDCLASS

# The inside() class functions are used to test whether a given mesh
# point lies inside the boundary subdomain.

#STARTCLASS
class XBoundary_C(df_m.SubDomain):
    """The XBoundary_C class is a specialized SubDomain
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

#class XBoundary_C(df_m.SubDomain):ENDCLASS

#STARTCLASS
class YBoundary_C(df_m.SubDomain):
    """The YBoundary_C class is a specialized SubDomain
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

#class YBoundary_C(df_m.SubDomain):ENDCLASS

#STARTCLASS
class SourceXrange(df_m.SubDomain):
    """The SourceXrange class is a specialized SubDomain to used to
       mark points that are inside this source region.

    """
    # Set the X range of the source region in the constructor
    def __init__(self, x_min_value, x_max_value, map_tol=1.0e-10):
        self.x_min_value=x_min_value
        self.x_max_value=x_max_value
        # Pass the map_tol argument on to the base class
        super(self.__class__, self).__init__(map_tol=1.0e-10)

    # Return true for points inside the subdomain
    def inside(self, x, on_boundary):
        tol = 1.0e-10
        return x[0] > self.x_min_value-tol and x[0] < self.x_max_value+tol

#class SourceXrange(df_m.SubDomain):ENDCLASS

# User exposes whatever mesh parameters are useful in __init__ and
# these can be set in __main__

#STARTCLASS
class UserMesh_C(Mesh_C):
    """UserMesh_C is derived from Mesh_C and creates a 1D, 2D, or 3D Cartesian mesh.

       The units are MKS by default (i.e., if no conversion factor is applied), but
       any units can be used provided the conversion to MKS is available in the
       UserUnits_M.py module.

       Methods: 
           compute_cell_index() : Inherited from Mesh_C
           create_mesh() : Defined in this subclass
    """

# Select the unit system to be used for input parameters.
#    Convert = U_M.MyPlasmaUnits_C

    # Constructor
    def __init__(self, mesh_input, compute_dictionaries=False, compute_tree=False, plot_flag=False, plot_title=None):
        """
            The class UserMesh_C contains these attributes:
                1. A mesh.
                2. A MeshFunction that marks mesh boundaries for boundary conditions.
        """

        self.coordinate_system = 'Cartesian'

        self.create_mesh(mesh_input, plot_flag=plot_flag, plot_title=plot_title)

        # Call the parent constructor to complete setting class variables.
        super(UserMesh_C, self).__init__(mesh_file=None, compute_dictionaries=compute_dictionaries, compute_tree=compute_tree, plot_flag=plot_flag, plot_title=plot_title)

        self.field_boundary_dict = mesh_input.field_boundary_dict
        self.particle_boundary_dict = mesh_input.particle_boundary_dict

        return

    # Methods

#class UserMesh_C(Mesh_C):
    def create_mesh(self, mesh_input, plot_flag=False, plot_title=None):
        """
           Create a mesh according to the user's specifications.
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

        umi = mesh_input

#        print 'umi.cells_on_side = ', umi.cells_on_side

        xmin = umi.pmin.x()
        xmax = umi.pmax.x()
        ymin = umi.pmin.y()
        ymax = umi.pmax.y()
        zmin = umi.pmin.z()
        zmax = umi.pmax.z()

        # Options: 'left, 'right', 'left/right', 'crossed'
        diagonal = umi.diagonal

        # Boundary conditions for fields and particles
        fieldBoundaryDict = umi.field_boundary_dict
        particleBoundaryDict = umi.particle_boundary_dict
        particleSourceDict = umi.particle_source_dict

        plotTitle = plot_title

#        if umi.pmin.y() == umi.pmax.y() and umi.pmin.z() == umi.pmax.z():
        if ymin == ymax and zmin == zmax:
            # 1D mesh
            (nx,) = umi.cells_on_side # Need the comma to indicate a tuple
            mesh_df = df_m.IntervalMesh(nx, xmin, xmax)

            if plot_title is None: plotTitle = "X-mesh"
            # Name the mesh file that will be written below
            meshFileName = "X-mesh.pvd"

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
                fieldBoundaryMarker = df_m.MeshFunction('size_t', mesh_df, mesh_df.topology().dim()-1)

    #        print 'dim =', mesh_df.topology().dim()

                # Initialize all the facets with a default value of 0.
                # Then overwrite boundary facets that have Dirichlet BCs.
                fieldBoundaryMarker.set_all(0)

                # Create mesh subdomains Gamma_nnn that are the boundaries
                # (Could also do this later when the mesh has been stretched)
                Gamma_xmin = XBoundary_C(xmin) # Gamma_xmin is the lower X boundary
                Gamma_xmax = XBoundary_C(xmax) # Gamma_xmax is the upper X boundary

                # Assign integer IDs to the boundaries with Dirichlet
                # conditions.  Create a FacetFunction to store these values.
                xmin_indx = fieldBoundaryDict['xmin']
                xmax_indx = fieldBoundaryDict['xmax']
                Gamma_xmin.mark(fieldBoundaryMarker, xmin_indx) # Mark the lower
                                                                # X boundary
                                                                # with xmin_indx
                Gamma_xmax.mark(fieldBoundaryMarker, xmax_indx) # Mark the upper
                                                                # X boundary
                                                                # with xmax_indx
            else:
                fieldBoundaryMarker = None
            #   END:if fieldBoundaryDict is not None:

            if particleBoundaryDict is not None:
                # Mark the facets for particle boundary conditions
                # Create a MeshFunction object defined on the facets of this
                # mesh.  The function has a default (size_t) value of 0,
                # meaning 'no action needed'

                # A boundary has dimension 1 less than the domain:
                particleBoundaryMarker = df_m.MeshFunction('size_t', mesh_df, mesh_df.topology().dim()-1)
                # Initialize all mesh facets with a default value of 0.
                particleBoundaryMarker.set_all(0)
            else:
                particleBoundaryMarker = None
            #   END:if particleBoundaryDict is not None:

            if particleSourceDict is not None:
                # Mark the cells for particle source regions.
                # Create a MeshFunction object defined on the mesh cells
                # The function has a default (size_t) value of 0,
                # meaning 'no particle creation here'

                particleSourceMarker = df_m.CellFunction('size_t', mesh_df)
                particleSourceMarker.set_all(0)
                # Create mesh subdomains where particles are generated

                # There's one line below for each entry in the dictionary
                sourceX1 = SourceXrange(sourceX1min, sourceX1max)
                sourceX2 = SourceXrange(sourceX2min, sourceX2max)

                particleSourceDictInv = {v: k for k, v in particleSourceDict.items()}
                # Get the numerical value that identifies this source
                sourceX1_indx = int(particleSourceDictInv['sourceX1'])
                sourceX2_indx = int(particleSourceDictInv['sourceX2'])
                # Mark the cells in the source with this value
                sourceX1.mark(particleSourceMarker, sourceX1_indx)
                sourceX2.mark(particleSourceMarker, sourceX2_indx)
            else:
                particleSourceMarker = None
            #   END:if particleSourceDict is not None:

#        elif umi.pmin.z() == umi.pmax.z():
        elif zmin == zmax:
            # 2D mesh
            (nx, ny) = umi.cells_on_side

# v > 1.5:
#            if df_m.DOLFIN_VERSION_STRING > '1.5.0':
            if df_m.__version__ > '1.5.0':
                if diagonal is not None:
                    mesh_df = df_m.RectangleMesh(umi.pmin, umi.pmax, nx, ny, diagonal)
                else:
                    mesh_df = df_m.RectangleMesh(umi.pmin, umi.pmax, nx, ny)
            else:
# v = 1.5:
                mesh_df = df_m.RectangleMesh(umi.pmin[0], umi.pmin[1], umi.pmax[0], umi.pmax[1], nx, ny)

            if plot_title is None: plotTitle = "XY-mesh"
            # Name the mesh file that will be written below
            meshFileName = "XY-mesh.pvd"

            if fieldBoundaryDict is not None:
                error_msg = "Error in %s" % fncName
                sys.exit(error_msg)
            else:
                fieldBoundaryMarker = None

            if particleBoundaryDict is not None:
                # A boundary has dimension 1 less than the domain:
                particleBoundaryMarker = df_m.MeshFunction('size_t', mesh_df, mesh_df.topology().dim()-1)
                # Initialize all mesh facets with a default value of 0.
                particleBoundaryMarker.set_all(0)

                # Create mesh subdomains to apply boundary-conditions
                # There's one line below for each entry in the dictionary
                Gamma_xmin = XBoundary_C(xmin)
                Gamma_xmax = XBoundary_C(xmax)
                Gamma_ymin = YBoundary_C(ymin)
                Gamma_ymax = YBoundary_C(ymax)

#                print "particleBoundaryDict =", particleBoundaryDict

                # Mark these subdomains (boundaries) with non-zero size_t integers

                xmin_indx = particleBoundaryDict['xmin']
                xmax_indx = particleBoundaryDict['xmax']
                ymin_indx = particleBoundaryDict['ymin']
                ymax_indx = particleBoundaryDict['ymax']
                Gamma_xmin.mark(particleBoundaryMarker, xmin_indx)
                Gamma_xmax.mark(particleBoundaryMarker, xmax_indx)
                Gamma_ymin.mark(particleBoundaryMarker, ymin_indx)
                Gamma_ymax.mark(particleBoundaryMarker, ymax_indx)
            else:
                particleBoundaryMarker = None
#           END:if particleBoundaryDict is not None:

            if particleSourceDict is not None:
                error_msg = "Error in %s" % fncName
                sys.exit(error_msg)
            else:
                particleSourceMarker = None

#       END:elif zmin == zmax:

        else:
            # 3D mesh
            (nx, ny, nz) = umi.cells_on_side
# v > 1.5
            if df_m.__version__ > '1.5.0':
                mesh_df = df_m.BoxMesh(umi.pmin, umi.pmax, nx, ny, nz)
            else:
# v = 1.5:
                mesh_df = df_m.BoxMesh(umi.pmin[0], umi.pmin[1], umi.pmin[2], umi.pmax[0], umi.pmax[1], umi.pmax[2], nx, ny, nz)

            if plot_title is None: plotTitle = "XYZ-mesh"
            # Name the mesh file that will be written below
            meshFileName = "XYZ-mesh.pvd"

            if fieldBoundaryDict is not None:
                error_msg = "Error in %s" % fncName
                sys.exit(error_msg)
            else:
                fieldBoundaryMarker = None

            if particleBoundaryDict is not None:
            # A boundary has dimension 1 less than the domain:
                particleBoundaryMarker = df_m.MeshFunction('size_t', mesh_df, mesh_df.topology().dim()-1)
            # Initialize all mesh facets with a default value of 0.
                particleBoundaryMarker.set_all(0)
            else:
                particleBoundaryMarker = None

            if particleSourceDict is not None:
                error_msg = "Error in %s" % fncName
                sys.exit(error_msg)
            else:
                particleSourceMarker = None

#       END:if ymin == ymax and zmin == zmax:

        # Make a plot of the mesh, with non-zero values showing marked
        # boundaries
        if (plot_flag):
            df_m.plot(mesh_df, title=plotTitle, axes=True)
#            df_m.plot(fieldBoundaryMarker, title='field boundary marks', axes=True)
            mplot_m.show()
#            yesno = raw_input("Just called show() in create_mesh")

#            df_m.plot(particleBoundaryMarker, title='particle boundary marks', axes=True)
#            df_m.interactive()

        # Write the mesh to a VTK file for plotting
        meshFile = df_m.File(meshFileName)
        meshFile << mesh_df

        # Save the class attributes
        self.mesh = mesh_df
        self.field_boundary_marker = fieldBoundaryMarker
        self.particle_boundary_marker = particleBoundaryMarker
        self.particle_source_marker = particleSourceMarker

        return
#    def create_mesh(self, mesh_input):ENDDEF

#class UserMesh_C(Mesh_C):ENDCLASS

#
# Solve the equations for the fields
#
class UserPoissonSolve1D_C(PoissonSolve_C):
    """UserPoissonSolve1D_C solves the 1D Poisson equation in Cartesian coordinates to
       obtain fields from sources.  

       It can be edited by the user to specify the field solvers and solver tuning
       parameters, if any.  The units are MKS by default (i.e., if no conversion
       factor is applied), but any units can be used provided the conversion to MKS
       is available in the UserUnits_M.py module.
    """

    def __init__(self, phi_F, linear_solver, preconditioner, field_boundary_marker, phi_BCs, neg_electric_field=None):
        """mesh argument is only needed if, e.g., using a SpatialCoordinate in the equations.
        """

        # Select the unit system to be used for input parameters.
        constants = U_M.MyPlasmaUnits_C

        # Call the PoissonSolve_C base class constructor for
        # non-problem-specific initialization.
        super(self.__class__, self).__init__(phi_F)

        ### Set class and local variables from the arguments ###

       # Field-solver parameters
        self.solver_parameters = {}
        if linear_solver is not None:
            self.solver_parameters['linear_solver'] = linear_solver
        if preconditioner is not None:
            self.solver_parameters['preconditioner'] = preconditioner

        self.neg_electric_field = neg_electric_field

        # Extract the Dirichlet boundary conditions
        (xmin_indx, phi_xmin) = phi_BCs['xmin']
        (xmax_indx, phi_xmax) = phi_BCs['xmax']

        # Create functions from the boundary values        
        u_xmin = df_m.Constant(phi_xmin)
        u_xmax = df_m.Constant(phi_xmax)

        ### Create the Dirichlet boundary-condition list for the Poisson PDE.

        # The indices xmin_indx and xmax_indx identify the
        # facets in boundary_marker where boundary values are to be
        # set.
        #       args: DirichletBC(FunctionSpace, GenericFunction, MeshFunction, int, method="topological")

        self.bcs = [df_m.DirichletBC(self.V, u_xmin, field_boundary_marker, xmin_indx), df_m.DirichletBC(self.V, u_xmax, field_boundary_marker, xmax_indx)]

        ### Set up the variational problem ###

        ## Make the bilinear form 'a(w,v)' ##

        # The bilinear form 'a' is the LHS in the variational form of the
        # Laplace/Poisson eq. It is the 'inner product' of grad(w) times
        # grad(v), i.e., the integral of this product (times any coefficients)
        # over the domain. The spherical-coordinate form has a coefficient r**2
        # (for the volume-element needed to integrate over the domain). This
        # Cartesian-coordinate form doesn't.
        epsilon = df_m.Constant(constants.epsilon_0)
        self.a = epsilon*df_m.inner(df_m.nabla_grad(self.w), df_m.nabla_grad(self.v))*df_m.dx

        # Specify whether 'a(u,v)' has time-independent coefficients.  If so, the
        # matrix A only needs to be assembled once, avoiding redundant work.
        # By inspection of the expression above:
        self.pde_has_constant_coeffs = True

        # Carry out the initial assembly
        self.assemble_matrix()

        # Initialize the source vector 'b' to zero charge-density.
        self.assemble_source_expression(0.0)
        
        # Set the level of diagnostic output from the solver.
#        df_m.set_log_level(df_m.PROGRESS) # df.set_log_level(1) gives the most messages
        df_m.set_log_level(df_m.LogLevel.PROGRESS) # df.set_log_level(1) gives the most messages
#        df_m.set_log_level(16) # Same as PROGRESS level

# default LU is flakey: different answers on different calls: NO!: this was a heap problem of unitialized memory!
#        self.phi = None
#        self.negE = None

        return
#    def __init__(self, phi, linear_solver, preconditioner, field_boundary_marker, phi_BCs, charge_density=None, neg_electric_field=None):ENDDEF

#class UserPoissonSolve1D_C(PoissonSolve_C):ENDCLASS
