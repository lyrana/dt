# Make a mesh that DT can use

#__version__ = 0.1
#__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

"""UserMesh defines the mesh.
"""

import sys
import os
import numpy as np_M
import math
# !!! Direct invocation of dolfin. OK because UserMesh_C is a
# sub-class of Mesh_C !!!
import dolfin as df_M

from Dolfin_Module import Mesh_C

import UserUnits_Module as U_M

### Define boundary geometry

# The inside() class functions are used to test whether a given mesh
# point lies inside the boundary subdomain.

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

### Define the geometry of the particle sources

class SourceXrange(df_M.SubDomain):
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
        return True if (x[0] > self.x_min_value-tol and x[0] < self.x_max_value+tol) else False

# User exposes whatever mesh parameters are useful in __init__ and
# these can be set in __main__

class UserMesh_C(Mesh_C):
    """UserMesh_C is derived from Mesh_C.  It is to be edited by the user to specify the simulation
       mesh.  The units are MKS by default (i.e., if no conversion
       factor is applied), but any units can be used provided the
       conversion to MKS is available in the UserUnits_M.py module.

       This class creates a 2D or 3D rectangular mesh.  It was written
       to test particle cell indexing, so no boundary-conditions are
       set.

       Methods: 
           compute_cell_index() : Inherited from Mesh_C
           create_mesh() : Defined in this subclass
    """

# Select the unit system to be used for input parameters.
    Convert = U_M.MyPlasmaUnits_C

    # Constructor
    def __init__(self, meshInputCI=None, compute_dictionaries=False, compute_tree=False, plot_flag=False):
        """
            The class UserMesh_C contains these attributes:
                1. A mesh.
                2. A MeshFunction that marks mesh boundaries for boundary conditions.
        """

        if meshInputCI is not None:
            self.create_mesh(meshInputCI)
        else:
            error_msg = "UserMesh_C: no mesh input given"
            sys.exit(error_msg)

        # Call the parent constructor to complete setting class variables.
        super(UserMesh_C, self).__init__(mesh_file=None, compute_dictionaries=compute_dictionaries, compute_tree=compute_tree, plot_flag=plot_flag)

        self.field_boundary_dict = meshInputCI.field_boundary_dict
        self.particle_boundary_dict = meshInputCI.particle_boundary_dict

        return

    # Methods

#class UserMesh_C(Mesh_C):
    def create_mesh(self, meshInputCI):
        """
           Create a mesh according to the user's specifications.
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

        miCI = meshInputCI

        xmin = miCI.pmin.x()
        xmax = miCI.pmax.x()
        ymin = miCI.pmin.y()
        ymax = miCI.pmax.y()
        zmin = miCI.pmin.z()
        zmax = miCI.pmax.z()

        # Boundary conditions for fields and particles
        fieldBoundaryDict = miCI.field_boundary_dict
        particleBoundaryDict = miCI.particle_boundary_dict
        particleSourceDict = miCI.particle_source_dict

#        if miCI.pmin.y() == miCI.pmax.y() and miCI.pmin.z() == miCI.pmax.z():
        if ymin == ymax and zmin == zmax:
            # 1D mesh
            (nx) = miCI.cells_on_side
            mesh = df_M.IntervalMesh(nx, xmin, xmax)

            # Mark the facets for particle boundary conditions
            # Create a MeshFunction object defined on the facets of this
            # mesh.  The function has a default (size_t) value of 0,
            # meaning 'no action needed'

            # A boundary has dimension 1 less than the domain:
            particleBoundaryMarker = df_M.MeshFunction('size_t', mesh, mesh.topology().dim()-1)
            # Initialize all mesh facets with a default value of 0.
            particleBoundaryMarker.set_all(0)

            # Mark the cells for particle source regions.
            # Create a MeshFunction object defined on the mesh cells
            # The function has a default (size_t) value of 0,
            # meaning 'no particle creation here'

            particleSourceMarker = df_M.CellFunction('size_t', mesh)
            particleSourceMarker.set_all(0)

            if particleSourceDict is not None:
                # Create mesh subdomains where particles are generated

                # There's one line below for each entry in the dictionary
                sourceX1 = SourceXrange(sourceX1min, sourceX1max)
                sourceX2 = SourceXrange(sourceX2min, sourceX2max)

                particleSourceDictInv = {v: k for k, v in particleSourceDict.iteritems()}
                # Get the numerical value that identifies this source
                sourceX1_indx = int(particleSourceDictInv['sourceX1'])
                sourceX2_indx = int(particleSourceDictInv['sourceX2'])
                # Mark the cells in the source with this value
                sourceX1.mark(particleSourceMarker, sourceX1_indx)
                sourceX2.mark(particleSourceMarker, sourceX2_indx)

#        elif miCI.pmin.z() == miCI.pmax.z():
        elif zmin == zmax:
            # 2D mesh
            (nx, ny) = miCI.cells_on_side
# v > 1.5:
            if df_M.DOLFIN_VERSION_STRING > '1.5.0':
                mesh = df_M.RectangleMesh(miCI.pmin, miCI.pmax, nx, ny)
            else:
# v = 1.5:
                mesh = df_M.RectangleMesh(miCI.pmin[0], miCI.pmin[1], miCI.pmax[0], miCI.pmax[1], nx, ny)

            # A boundary has dimension 1 less than the domain:
            particleBoundaryMarker = df_M.MeshFunction('size_t', mesh, mesh.topology().dim()-1)
            # Initialize all mesh facets with a default value of 0.
            particleBoundaryMarker.set_all(0)

            if particleBoundaryDict is not None:
                # Create mesh subdomains to apply boundary-conditions
                # There's one line below for each entry in the dictionary
                Gamma_xmin = XBoundary(xmin)
                Gamma_xmax = XBoundary(xmax)
                Gamma_ymin = YBoundary(ymin)
                Gamma_ymax = YBoundary(ymax)

#                print "particleBoundaryDict =", particleBoundaryDict

                # Mark these subdomains (boundaries) with non-zero integers

                xmin_indx = particleBoundaryDict['xmin']
                xmax_indx = particleBoundaryDict['xmax']
                ymin_indx = particleBoundaryDict['ymin']
                ymax_indx = particleBoundaryDict['ymax']
                Gamma_xmin.mark(particleBoundaryMarker, xmin_indx)
                Gamma_xmax.mark(particleBoundaryMarker, xmax_indx)
                Gamma_ymin.mark(particleBoundaryMarker, ymin_indx)
                Gamma_ymax.mark(particleBoundaryMarker, ymax_indx)
#           END:if particleBoundaryDict is not None:
#       END:elif zmin == zmax:

        else:
            # 3D mesh
            (nx, ny, nz) = miCI.cells_on_side
# v > 1.5
            if df_M.DOLFIN_VERSION_STRING > '1.5.0':
                mesh = df_M.BoxMesh(miCI.pmin, miCI.pmax, nx, ny, nz)
            else:
# v = 1.5:
                mesh = df_M.BoxMesh(miCI.pmin[0], miCI.pmin[1], miCI.pmin[2], miCI.pmax[0], miCI.pmax[1], miCI.pmax[2], nx, ny, nz)
            # A boundary has dimension 1 less than the domain:
            particleBoundaryMarker = df_M.MeshFunction('size_t', mesh, mesh.topology().dim()-1)
            # Initialize all mesh facets with a default value of 0.
            particleBoundaryMarker.set_all(0)
#       END:if ymin == ymax and zmin == zmax:

        # Save the class attributes
        self.mesh = mesh
        self.particle_boundary_marker = particleBoundaryMarker

        return
#    def create_mesh(self, meshInputCI):ENDDEF

#class UserMesh_C(Mesh_C):ENDCLASS
