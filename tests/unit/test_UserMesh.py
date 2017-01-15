#!/usr/bin/env python

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import os
import math
import unittest

import dolfin as df_M

from DT_Module import DTmeshInput_C
from UserMesh_y_Fields_FE2D_Module import UserMesh_C

class TestUserMesh_y_Fields(unittest.TestCase):
    """Test the mesh class in UserMesh_y_Fields_FE2D_Module.py"""
    
    def setUp(self):
        # initializations for each test go here...

        miCI = DTmeshInput_C()

        # Make the mesh
        # radial
        miCI.rmin, miCI.rmax = 1.0, 5.0 # Mesh goes from rmin to rmax in radius
        miCI.nr = 10 # Number of divisions in r direction
        miCI.stretch = 1.3 # Stretch parameter

        # theta, starts at 0
        miCI.tmax = math.pi/2 # quarter-circle
        miCI.nt = 20  # Number of divisions in theta direction
        
        # The diagonal that makes the triangular mesh
        # Options: 'left, 'right', 'left/right', 'crossed'
        miCI.diagonal = 'crossed'

        # Name the Dirichlet boundaries and assign integers to them.
        # These are the boundary-name -> int pairs used to mark mesh
        # facets:
        rmin_indx = 1
        rmax_indx = 2
        fieldBoundaryDict = {'rmin': rmin_indx,
                             'rmax': rmax_indx,
                             }

        miCI.field_boundary_dict = fieldBoundaryDict

        # Name the boundaries used to apply particle
        # boundary-conditions and assign integers to them.  These are
        # the boundary-name -> int pairs used to mark particle-mesh
        # facets:
        rmin_indx = 1
        rmax_indx = 2
        thmin_indx = 4
        thmax_indx = 8
        particleBoundaryDict = {'rmin': rmin_indx,
                                'rmax': rmax_indx,
                                'thmin': thmin_indx,
                                'thmax': thmax_indx,
                                }

        miCI.particle_boundary_dict = particleBoundaryDict

        self.miCI = miCI

        return

    def test_quarter_circle_plot_false(self):

        fncname = sys._getframe().f_code.co_name
        print '\ntest: ', fncname, '('+__file__+')'

        meshCI = UserMesh_C(meshInputCI=self.miCI, plotFlag=False)

#        df_M.plot(meshCI.mesh, title='cylindrical mesh', axes=True)
#        df_M.interactive()

#        yesno = raw_input("Looks OK [Y/n]?")
#        self.assertNotEqual(yesno, 'n', "Problem with mesh")
        
    def test_quarter_circle_plot_true(self):

        fncname = sys._getframe().f_code.co_name
        print '\ntest: ', fncname, '('+__file__+')'

        if os.environ.get('DISPLAY') is None:
            plotFlag=False
        else:
            plotFlag=True

        meshCI = UserMesh_C(meshInputCI=self.miCI, plotFlag=plotFlag)

#        df_M.plot(meshCI.mesh, title='cylindrical mesh', axes=True)
#        df_M.interactive()

#        yesno = raw_input("Looks OK [Y/n]?")
#        self.assertNotEqual(yesno, 'n', "Problem with mesh")

        # Write the mesh to a file:
        mesh_file = df_M.File('quarter_circle_mesh_crossed.xml') # Could use if-test on the value of 'diagonal'
        mesh_file << meshCI.mesh
        # Read back in:
        mesh_file >> meshCI.mesh

        # Write the boundary marker function to a file:
        field_boundary_marker_file = df_M.File('Fbcs_quarter_circle_mesh_crossed.xml') # Could use if-test on the value of 'diagonal'
        field_boundary_marker_file << meshCI.field_boundary_marker
        # Read back in:
        field_boundary_marker_file >> meshCI.field_boundary_marker

        particle_boundary_marker_file = df_M.File('Pbcs_quarter_circle_mesh_crossed.xml') # Could use if-test on the value of 'diagonal'
        particle_boundary_marker_file << meshCI.particle_boundary_marker
        # Read back in:
        particle_boundary_marker_file >> meshCI.particle_boundary_marker

if __name__ == '__main__':
    unittest.main()
