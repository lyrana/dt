#!/usr/bin/env python3

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import os
import math
import unittest

import dolfin as df_m


class TestUserMesh(unittest.TestCase):
    """Test the mesh class in UserMesh_y_Fields modules"""
    
#class TestUserMesh(unittest.TestCase):
    def setUp(self):
        # initializations for each test go here...


        # 1D spherical-coordinate mesh
        from UserMesh_y_Fields_Spherical1D_Module import UserMeshInput1DS_C, UserMesh1DS_C

        umi1DS = UserMeshInput1DS_C()

        ## Make the mesh

        # radial, r = 1 to 5 (Changed in one of the tests to r = 0 to 4)
        umi1DS.rmin, umi1DS.rmax = 1.0, 5.0 # Mesh goes from rmin to rmax in radius
        umi1DS.nr = 10 # Number of divisions in r direction
        umi1DS.stretch = 1.3 # Stretch parameter

        # Name the Dirichlet boundaries and assign integers to them.
        # These are the boundary-name -> int pairs used to mark mesh
        # facets:
        rminIndx = 1
        rmaxIndx = 2
        fieldBoundaryDict = {'rmin': rminIndx,
                             'rmax': rmaxIndx,
                             }

        umi1DS.field_boundary_dict = fieldBoundaryDict

        # Name the boundaries used to apply particle
        # boundary-conditions and assign integers to them.  These are
        # the boundary-name -> int pairs used to mark particle-mesh
        # facets:
        rmin_indx = 1
        rmax_indx = 2
        particleBoundaryDict = {'rmin': rmin_indx,
                                'rmax': rmax_indx,
                                }

        umi1DS.particle_boundary_dict = particleBoundaryDict

        self.umi1DS = umi1DS

        # 2D quarter-circle mesh
        from UserMesh_y_Fields_FE2D_Module import UserMeshInput2DCirc_C

        umi2DCirc = UserMeshInput2DCirc_C()

        # Make the mesh
        # radial
        umi2DCirc.rmin, umi2DCirc.rmax = 1.0, 5.0 # Mesh goes from rmin to rmax in radius
        umi2DCirc.nr = 10 # Number of divisions in r direction
        umi2DCirc.stretch = 1.3 # Stretch parameter

        # theta, starts at 0
        umi2DCirc.tmax = math.pi/2 # quarter-circle
        umi2DCirc.nt = 20  # Number of divisions in theta direction
        
        # The diagonal that makes the triangular mesh
        # Options: 'left, 'right', 'left/right', 'crossed'
        umi2DCirc.diagonal = 'crossed'

        # Name the Dirichlet boundaries and assign integers to them.
        # These are the boundary-name -> int pairs used to mark mesh
        # facets:
        rmin_indx = 1
        rmax_indx = 2
        fieldBoundaryDict = {'rmin': rmin_indx,
                             'rmax': rmax_indx,
                             }

        umi2DCirc.field_boundary_dict = fieldBoundaryDict

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

        umi2DCirc.particle_boundary_dict = particleBoundaryDict

        self.umi2DCirc = umi2DCirc

        return

#class TestUserMesh(unittest.TestCase):
    def test_1D_spherical_mesh(self):
        """At the moment, this just writes mesh files for other tests.

           The mesh starts at r = 1.
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName, '('+__file__+')')

        from UserMesh_y_Fields_Spherical1D_Module import UserMesh1DS_C

        if os.environ.get('DISPLAY') is None:
            plotFlag=False
        else:
            plotFlag=True

        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name
        mesh_M = UserMesh1DS_C(self.umi1DS, compute_tree=False, plot_flag=plotFlag, plot_title=plotTitle)

        # Write the mesh to a file:
        mesh_file = df_m.File('mesh_1D_radial.xml')
        mesh_file << mesh_M.mesh

        # Write the boundary marker function to a file:
        # Field boundaries
        field_boundary_marker_file = df_m.File('mesh_1D_radial_Fbcs.xml')
        field_boundary_marker_file << mesh_M.field_boundary_marker

        # Particle boundaries
        particle_boundary_marker_file = df_m.File('mesh_1D_radial_Pbcs.xml')
        particle_boundary_marker_file << mesh_M.particle_boundary_marker

        return
#    def test_1D_spherical_mesh(self):ENDDEF

#class TestUserMesh(unittest.TestCase):
    def test_1D_spherical_mesh_r0(self):
        """At the moment, this just writes mesh files for other tests.

           The mesh starts at r = 0.
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName, '('+__file__+')')

        from UserMesh_y_Fields_Spherical1D_Module import UserMesh1DS_C

        if os.environ.get('DISPLAY') is None:
            plotFlag=False
        else:
            plotFlag=True


        self.umi1DS.rmin, self.umi1DS.rmax = 0.0, 4.0 # Mesh goes from rmin to rmax in radius
        
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name
        mesh_M = UserMesh1DS_C(self.umi1DS, compute_tree=False, plot_flag=plotFlag, plot_title=plotTitle)

        # Write the mesh to a file:
        mesh_file = df_m.File('mesh_1D_radial_r0.xml')
        mesh_file << mesh_M.mesh

        # Write the boundary marker function to a file:
        # Field boundaries
        field_boundary_marker_file = df_m.File('mesh_1D_radial_r0_Fbcs.xml')
        field_boundary_marker_file << mesh_M.field_boundary_marker

        # Particle boundaries
        particle_boundary_marker_file = df_m.File('mesh_1D_radial_r0_Pbcs.xml')
        particle_boundary_marker_file << mesh_M.particle_boundary_marker

        return
#    def test_1D_spherical_mesh_r0(self):ENDDEF

#class TestUserMesh(unittest.TestCase):
    def test_quarter_circle_plot_false(self):

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName, '('+__file__+')')

        from UserMesh_y_Fields_FE2D_Module import UserMesh2DCirc_C
        mesh_M = UserMesh2DCirc_C(self.umi2DCirc, plot_flag=False)

        return
#    def test_quarter_circle_plot_false(self):ENDDEF

#class TestUserMesh(unittest.TestCase):
    def test_quarter_circle_plot_true(self):

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName, '('+__file__+')')

        if os.environ.get('DISPLAY') is None:
            plotFlag=False
        else:
            plotFlag=True

        from UserMesh_y_Fields_FE2D_Module import UserMesh2DCirc_C
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name + ": mesh"
        mesh_M = UserMesh2DCirc_C(self.umi2DCirc, plot_flag=plotFlag, plot_title=plotTitle)

#        df_m.plot(mesh_M.mesh, title='cylindrical mesh', axes=True)
#        df_m.interactive()

#        yesno = raw_input("Looks OK [Y/n]?")
#        self.assertNotEqual(yesno, 'n', "Problem with mesh")

        # Write the mesh to a file:
        mesh_file = df_m.File('mesh_quarter_circle_crossed.xml') # Could use if-test on the value of 'diagonal'
        mesh_file << mesh_M.mesh
        # Read back in:
        mesh_file >> mesh_M.mesh

        # Write the boundary marker function to a file:
        field_boundary_marker_file = df_m.File('mesh_quarter_circle_crossed_Fbcs.xml') # Could use if-test on the value of 'diagonal'
        field_boundary_marker_file << mesh_M.field_boundary_marker
        # Read back in:
        field_boundary_marker_file >> mesh_M.field_boundary_marker

        particle_boundary_marker_file = df_m.File('mesh_quarter_circle_crossed_Pbcs.xml') # Could use if-test on the value of 'diagonal'
        particle_boundary_marker_file << mesh_M.particle_boundary_marker
        # Read back in:
        particle_boundary_marker_file >> mesh_M.particle_boundary_marker

        return
#    def test_quarter_circle_plot_true(self):ENDDEF

#class TestUserMesh(unittest.TestCase):ENDCLASS

if __name__ == '__main__':
    unittest.main()
