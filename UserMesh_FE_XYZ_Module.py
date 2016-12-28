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
           compute_cell_index() : inherited from Mesh_C
           create_mesh() : local
    """

# Select the unit system to be used for input parameters.
    Convert = U_M.MyPlasmaUnits_C

    # Constructor
    def __init__(self, meshInputCI=None, computeDictionaries=False, computeTree=False, plotFlag=False):
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
        super(UserMesh_C, self).__init__(meshFile=None, computeDictionaries=computeDictionaries, computeTree=computeTree, plotFlag=plotFlag)

        return

    # Methods
    def create_mesh(self, miCI):
        """
           Create a mesh according to the user's specifications.
        """
        fncname = sys._getframe().f_code.co_name

        if miCI.pmin.y() == miCI.pmax.y() and miCI.pmin.z() == miCI.pmax.z():
            (nx) = miCI.cells_on_side
            self.mesh = df_M.IntervalMesh(nx, miCI.pmin.x(), miCI.pmax.x())
        elif miCI.pmin.z() == miCI.pmax.z():
            (nx, ny) = miCI.cells_on_side
# v > 1.5:
            if df_M.DOLFIN_VERSION_STRING > "1.5.0":
                self.mesh = df_M.RectangleMesh(miCI.pmin, miCI.pmax, nx, ny)
            else:
# v = 1.5:
                self.mesh = df_M.RectangleMesh(miCI.pmin[0], miCI.pmin[1], miCI.pmax[0], miCI.pmax[1], nx, ny)
        else:
            (nx, ny, nz) = miCI.cells_on_side
# v > 1.5
            if df_M.DOLFIN_VERSION_STRING > "1.5.0":
                self.mesh = df_M.BoxMesh(miCI.pmin, miCI.pmax, nx, ny, nz)
            else:
# v = 1.5:
                self.mesh = df_M.BoxMesh(miCI.pmin[0], miCI.pmin[1], miCI.pmin[2], miCI.pmax[0], miCI.pmax[1], miCI.pmax[2], nx, ny, nz)

        return
#    def create_mesh(self, miCI):ENDDEF

#class UserMesh_C(Mesh_C):ENDCLASS
