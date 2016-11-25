#!/usr/bin/env python

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

"""Get an error message from PETSc when no PETSc functions have been
   called explicitly.  Importing dolfin may be the reason.  Perhaps
   PETSc gets initialized.  Only seems to happen in IDLE.

   Coding comes from unit/test_UserMesh.py
"""
import math

import dolfin as df_M

from UserMesh_rtheta2D import Mesh_C

# radial
rmin, rmax = 1.0, 5.0 # Mesh goes from rmin to rmax in radius
nr = 10  # Number of divisions in r direction
s = 1.3 # Stretch parameter

# theta, starts at 0
thetamax = math.pi/2 # quarter-circle
nt = 20  # Number of divisions in theta direction

#meshObj = Mesh_C(rmin, rmax, s, nr, thetamax, nt, plotFlag=False)
meshObj = Mesh_C(rmin, rmax, s, nr, thetamax, nt, plotFlag=True)
df_M.plot(meshObj.mesh, title='cylindrical mesh', axes=True)
df_M.interactive()
yesno = raw_input("Looks OK [Y/n]?")

