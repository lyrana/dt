#!/usr/bin/env python3

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import os
import math
import unittest

import dolfin as df_m
import matplotlib.pyplot as mplot_m

from UserMesh_y_Fields_FE_XYZ_Module import *

class TestMeshFunctions(unittest.TestCase):
    """Test algorithms that create lookup tables for meshes.
    
     """
    
    def setUp(self):
        # initializations for each test go here...

        self.plotMesh = False
        self.plotResults = False

        # Turn plots off if there's no display.
        if os.environ.get('DISPLAY') is None:
            self.plotMesh = False
            self.plotResults = False

        return


    def test_1_cell_neighbors_list(self):
        """ For each cell in a mesh, make a list of its neighbor cells.

            Check the list against a list created by visual inspection of the mesh.

            This example is from:
              http://fenicsproject.org/qa/3843/cell-neighbours
        """

        mesh = df_m.UnitSquareMesh(2, 2)
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name + ": mesh"

        # Init facet-cell connectivity
        tdim = mesh.topology().dim()
        mesh.init(tdim - 1, tdim)

        expected_neighbors = {0: [3, 1], 1: [4, 0], 2: [3], 3: [6, 2, 0], 4: [7, 5, 1], 5: [4], 6: [7, 3], 7: [6, 4]}

        # For every cell, build a list of cells that are connected to its facets
        # but are not the iterated cell.

        # Create a dictionary, where the keys are the cell indices, and the values are a list of the cell indices of touching cells.
        # Double loop: outer one over cells in the mesh, inner one over facets of the cell.
        # cells(mesh) are the cells in the mesh. 
        # cell(index) is the local index of the cell.
        # facets(cell) are the facets of a cell.
        # facet.entities(tdim) are the entities of dimension 2 touching that facet
        cell_neighbors = {cell.index(): sum(([ci for ci in facet.entities(tdim) if ci != cell.index()] for facet in df_m.facets(cell)), []) for cell in df_m.cells(mesh)}

#        print cell vertex coordinates
#        for cell in df_m.cells(mesh):
#            print("cell index", cell.index(), "vertices", cell.get_vertex_coordinates())

        if self.plotMesh is True:
            df_m.plot(mesh, title=plotTitle)
            mplot_m.show()
        # Write the mesh to a VTK file
        meshFileName = "UnitSquareMesh_2by2.pvd"
        meshFile = df_m.File(meshFileName)
        meshFile << mesh            

        self.assertEqual(cell_neighbors, expected_neighbors, msg="Cell neighbor list is not correct")

        return


    def test_2_cell_facet_normals(self):
        """ Test the function compute_cell_facet_normals_dict().

            Check the list of normals against a manual computation of the normals.

            Note: This uses a different mesh from the previous test. File
            XY-8-tris-left-diag.png has a picture of the mesh.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName, '('+__file__+')')
        
        # Describe a 2D mesh from (-10,-10) to (10,10) with 2 cells on a side.
        # The mesh is triangular, so there's a total of 4x2 = 8 cells.
        umi2d_I = UserMeshInput_C()
        umi2d_I.pmin = df_m.Point(-10.0, -10.0)
        umi2d_I.pmax = df_m.Point(10.0, 10.0)
        umi2d_I.cells_on_side = (2, 2)
        umi2d_I.diagonal = 'left'

        # UserMesh_FE_XYZ_Module can make the mesh from the above input.
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name + ": XY mesh"
        mesh2d_M = UserMesh_C(umi2d_I, plot_flag=self.plotMesh, plot_title=plotTitle)
        
        # Compute the facet-normals:
        mesh2d_M.compute_cell_facet_normals_dict() # Unit vectors normal to cell facets        
        # print("cell_facet_normals_dict: ", mesh2d_M.cell_facet_normals_dict)
        
        cell_0_expected_facet_normals = [[0.70710678, 0.70710678], # opposite vertex 0
                                         [-1., 0.], # opposite vertex 1
                                         [ 0., -1.] # opposite vertex 3
                                        ]


        decimalPlaces = 6 # Accuracy of test
        # Check the computed values of the normal vector for the 3 facets of cell 0:
        for facet_id in range(3):
            for comp in range(2):
                self.assertAlmostEqual(mesh2d_M.cell_facet_normals_dict[0][facet_id][comp], cell_0_expected_facet_normals[facet_id][comp], places=decimalPlaces, msg="Cell facet normals for cell 0 are not correct")

        return
    
if __name__ == '__main__':
    unittest.main()
