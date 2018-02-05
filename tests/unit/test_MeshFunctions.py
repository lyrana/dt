#!/usr/bin/env python

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import os
import math
import unittest

import dolfin as df_m
import matplotlib.pyplot as mplot_m

from UserMesh_y_Fields_FE2D_Module import *

class TestMeshFunctions(unittest.TestCase):
    """Test the mesh class in UserMesh_y_Fields_FE2D_Module.py"""
    
    def setUp(self):
        # initializations for each test go here...

        return


    def test_1_cell_neighbor_list(self):
        """Example from:
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
        cell_neighbors = {cell.index(): sum((filter(lambda ci: ci != cell.index(), facet.entities(tdim)) for facet in df_m.facets(cell)), []) for cell in df_m.cells(mesh)}

#        print cell_neighbors

        df_m.plot(mesh, title=plotTitle)
        mplot_m.show()

#        self.assertEqual(len(shared_items), len(expected_neighbors), msg="Cell neighbor list is not correct")
        self.assertEqual(cmp(cell_neighbors, expected_neighbors), 0, msg="Cell neighbor list is not correct")

        return

if __name__ == '__main__':
    unittest.main()
