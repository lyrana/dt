#!/usr/bin/env python3

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import os
import math
import unittest

import numpy as np_m

import dolfin as df_m
import matplotlib.pyplot as mplot_m

from UserMesh_y_Fields_FE_XYZ_Module import *

# Use the C++ functions in the mesh_entity_arrays_solib.so library
import mesh_entity_arrays_solib


class TestMeshFunctions(unittest.TestCase):
    """Test C++ algorithms that create lookup tables for meshes.

       The lookup tables are members of the MeshEntityArrays class:
         cell_neighbors_array (accessed using 
          
    
     """
    
    def setUp(self):
        # initializations for each test go here...

        return


    def test_1_cell_neighbors_list(self):
        """For each cell in a 2D mesh, make a list of its neighbor-cells.

            The mesh is a unit square divided into 4x2 triangles. Check the
            cell-neighbor list against a list created by visual inspection of the
            mesh: UnitSquareMesh_2by2.png.

            This example is from:
              http://fenicsproject.org/qa/3843/cell-neighbours

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName, '('+__file__+')')
        
        # Create a 2D square mesh
        mesh_df = df_m.UnitSquareMesh(2, 2)
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name + ": mesh"

        # Create the name of the specialized MeshEntityArrays class with the right
        # number of cell-facets
        tDim = mesh_df.topology().dim()
        nFacets = tDim + 1
        meaClass = "MeshEntityArrays_" + str(nFacets) + "_facets"
        meaCtor = getattr(mesh_entity_arrays_solib, meaClass)
        
        # Call the MEA constructor
        meaObj = meaCtor(mesh_df, compute_particle_mesh_maps=True)        

        # The list from test_MeshFunctions.py:
        # expected_neighbors = {0: [3, 1], 1: [4, 0], 2: [3], 3: [6, 2, 0], 4: [7, 5, 1], 5: [4], 6: [7, 3], 7: [6, 4]}
        # Include the NO_CELL neighbors that are outside the mesh:
        expected_neighbors = {0: [3,1,-1], 1: [4,0,-1], 2: [-1,3,-1], 3: [6,2,0], 4: [7,5,1], 5: [-1,4,-1], 6: [-1,7,3], 7: [-1,6,4]}        

        # Get the neighbors from the C++ array
        cell_neighbors = {icell: meaObj.get_cell_neighbors(icell) for icell in range(mesh_df.num_cells())}
        
        df_m.plot(mesh_df, title=plotTitle)
        mplot_m.show()

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
        plotFlag = True
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name + ": XY mesh"
        mesh2d_M = UserMesh_C(umi2d_I, plot_flag=plotFlag, plot_title=plotTitle)
        mesh_df = mesh2d_M.mesh
        
        # Create the name of the specialized MeshEntityArrays class with the right
        # number of cell-facets
        tDim = mesh_df.topology().dim()
        nFacets = tDim + 1
        meaClass = "MeshEntityArrays_" + str(nFacets) + "_facets"
        meaCtor = getattr(mesh_entity_arrays_solib, meaClass)

        # Call the MEA constructor
        # print("Calling the MEA constructor...")
        meaObj = meaCtor(mesh_df, compute_particle_mesh_maps=True)        
        
        # Obtain the list of facet-normals for cell 0:
        # print("Obtain the facet-normals for cell 0...")
        cell_0_facet_normals_list = meaObj.get_cell_facet_normals(0) # Unit vectors normal to cell facets

        # Reshape the list to get nFacet vectors, with 3 components each:
        cell_0_facet_normals = np_m.reshape(cell_0_facet_normals_list, (nFacets, 3))
#        print("cell_0_facet_normals: ", cell_0_facet_normals)
        
        cell_0_expected_facet_normals = [[0.70710678, 0.70710678, 0.], # opposite vertex 0
                                         [-1., 0., 0.], # opposite vertex 1
                                         [ 0.,-1., 0.]  # opposite vertex 3
                                        ]

        # Check the computed values:
        decimalPlaces = 6 # Accuracy of test

        for facet_id in range(nFacets):
            for comp in range(3):
                self.assertAlmostEqual(cell_0_facet_normals[facet_id][comp], cell_0_expected_facet_normals[facet_id][comp], places=decimalPlaces, msg="Cell facet normals for cell 0 are not correct")

        return
    
if __name__ == '__main__':
    unittest.main()
