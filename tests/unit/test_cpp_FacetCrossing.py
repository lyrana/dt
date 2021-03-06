#!/usr/bin/env python3

__version__ = 0.1
__author__ = 'Copyright (C) 2019 L. D. Hughes'
#__all__ = []

import sys
import os
import numpy as np_m
import importlib as im_m
import unittest

import dolfin as df_m

from Dolfin_Module import Mesh_C
from UserUnits_Module import MyPlasmaUnits_C
from Particle_Module import *

from UserMesh_y_Fields_FE_XYZ_Module import *

# Use the C++ functions in the dolfin_functions_solib.so library
import dolfin_functions_solib


class TestCppFacetCrossing(unittest.TestCase):
    """Test the C++ algorithm to find the cell facet crossed and desination cell due
       to a particle displacement.

       Functions tested: find_facet()
    """
    
    def setUp(self):

        # Initializations performed before each test go here...

        plotFlag = False
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name

        # Turn off plotting if there's no DISPLAY

        # if os.environ.get('DISPLAY') is None:
        #     plotFlag=False
        # else:
        #     plotFlag=True

        # Create 1D, 2D and 3D meshes that the particles can be tested against

        # 1D mesh input
        umi1D = UserMeshInput_C()
        pmin = -0.03
        pmax = 0.03
        cells_on_side = 4
        umi1D.pmin = df_m.Point(pmin)
        umi1D.pmax = df_m.Point(pmax)
        umi1D.cells_on_side = (cells_on_side,) # Need the comma to indicate a tuple and
                                               # not just redundant parens.
        self.mesh1D_dx = (pmax-pmin)/cells_on_side

        # Create 1D mesh
        self.mesh1D = UserMesh_C(umi1D, compute_dictionaries=False, compute_cpp_arrays=True, compute_tree=True, plot_flag=plotFlag, plot_title=plotTitle + ": 1D")
        self.mesh1D.compute_cell_dict()
        self.mesh1D.compute_cell_entity_indices_dict('facet')        
        self.mesh1D.compute_cell_facet_normals_dict() # Unit vectors normal to cell facets

        # 2D mesh input
        umi2D = UserMeshInput_C()
        pmin = -0.03
        pmax = 0.03
        cells_on_side = 4
        umi2D.pmin = df_m.Point(pmin, pmin)
        umi2D.pmax = df_m.Point(pmax, pmax)
        umi2D.cells_on_side = (cells_on_side, cells_on_side)
#        umi2D.diagonal = 'crossed'
        self.mesh2D_dx = (pmax-pmin)/cells_on_side

        # Create 2D mesh
        self.mesh2D = UserMesh_C(umi2D, compute_dictionaries=False, compute_cpp_arrays=True, compute_tree=True, plot_flag=plotFlag, plot_title=plotTitle + ": 2D")
        self.mesh2D.compute_cell_dict()
        self.mesh2D.compute_cell_entity_indices_dict('facet')        
        self.mesh2D.compute_cell_facet_normals_dict() # Unit vectors normal to cell facets
        self.mesh2D.compute_cell_neighbors_dict()

        # 3D mesh input
        umi3D = UserMeshInput_C()
        umi3D.pmin = df_m.Point(-0.03, -0.03, -0.03)
        umi3D.pmax = df_m.Point(0.03, 0.03, 0.03)
        umi3D.cells_on_side = (4, 4, 4)
        # Create 3D mesh
        self.mesh3D = UserMesh_C(umi3D, compute_tree=True, plot_flag=plotFlag, plot_title=plotTitle + ": 3D")
        self.mesh3D.compute_cell_dict()
        # self.mesh3D.compute_cell_entity_indices_dict('vertex')
        self.mesh3D.compute_cell_entity_indices_dict('facet')
        self.mesh3D.compute_cell_facet_normals_dict() # Unit vectors normal to cell facets
        self.mesh3D.compute_cell_neighbors_dict()

        self.mesh3D.compute_cpp_arrays()

        # pmesh is the owner of the compute_index function?

        return
#    Setup():ENDDEF


    def test_1_facet_normals(self):
        """Test the facet-normal dictionary.
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest:', fncName, '('+__file__+')')

        # 1D mesh

        # To list the vertices of a facet, create facet-vertex
        # connectivity (both are dim 0 in 1D):
#        self.mesh1D.mesh.init(0,0)

        # Need following only if later code is uncommented:
        # vertex_coords = self.mesh1D.mesh.coordinates()

        nFacets = self.mesh1D.mesh.topology().dim() + 1
        for cell in df_m.cells(self.mesh1D.mesh):
            facetNormalsList = self.mesh1D.mea_object.get_cell_facet_normals(cell.index())
            # Reshape the list to get nFacet vectors, with 3 components each:
            facetNormalVectors = np_m.reshape(facetNormalsList, (nFacets, 3))

#            print("cell", cell.index(), "from MEA, facetNormalVectors", facetNormalVectors)

            mesh_dim = self.mesh1D.gdim
            for fi in range(mesh_dim+1):
                self.assertAlmostEqual(cell.normal(fi).x(), facetNormalVectors[fi,0], msg = "1D cell facet normal x component is not correct")

#            print "cell", cell.index(), "normal0", cell.normal(0).x(), "normal1", cell.normal(1).x()

        # 2D mesh
        nFacets = self.mesh2D.mesh.topology().dim() + 1        
        for cell in df_m.cells(self.mesh2D.mesh):
            facetNormalsList = self.mesh2D.mea_object.get_cell_facet_normals(cell.index())
            # Reshape the list to get nFacet vectors, with 3 components each:
            facetNormalVectors = np_m.reshape(facetNormalsList, (nFacets, 3))

            # Print values for a few cells
            for i in [0, 1, 6, 9,]:
                if cell.index() == i:
                    # use the "local index" 0, 1, 2 for the triangle facets
#                    print "cell", cell.index(), "normal0", cell.normal(0).x(), cell.normal(0).y()
#                    print "cell", cell.index(), "normal1", cell.normal(1).x(), cell.normal(1).y()
#                    print "cell", cell.index(), "normal2", cell.normal(2).x(), cell.normal(2).y()
                    pass

            mesh_dim = self.mesh2D.gdim
            for fi in range(mesh_dim+1):
                self.assertAlmostEqual(cell.normal(fi).x(), facetNormalVectors[fi,0], msg = "2D cell facet normal x component is not correct")
                self.assertAlmostEqual(cell.normal(fi).y(), facetNormalVectors[fi,1], msg = "2D cell facet normal y component is not correct")

        # 3D mesh
        nFacets = self.mesh3D.mesh.topology().dim() + 1        
        for cell in df_m.cells(self.mesh3D.mesh):
            facetNormalsList = self.mesh3D.mea_object.get_cell_facet_normals(cell.index())
            # Reshape the list to get nFacet vectors, with 3 components each:
            facetNormalVectors = np_m.reshape(facetNormalsList, (nFacets, 3))
            # Check all the normals for all cells
            mesh_dim = self.mesh3D.gdim
            for fi in range(mesh_dim+1):
                self.assertAlmostEqual(cell.normal(fi).x(), facetNormalVectors[fi,0], msg = "3D cell facet normal x component is not correct")
                self.assertAlmostEqual(cell.normal(fi).y(), facetNormalVectors[fi,1], msg = "3D cell facet normal y component is not correct")
                self.assertAlmostEqual(cell.normal(fi).z(), facetNormalVectors[fi,2], msg = "3D cell facet normal z component is not correct")

        return
#    def test_1_facet_normals(self):ENDDEF


    def test_2_facet_crossing(self):
        """ Test calculation of the cell facet crossed by a move vector.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest:', fncName, '('+__file__+')')

        #
        # 1D: Initial point = center of cell, move = cell length
        #
        numberOfFacetCrossings = 0
        for cell in df_m.cells(self.mesh1D.mesh):
            neighborCellIndices = self.mesh1D.mea_object.get_cell_neighbors(cell.index())

            # Create initial point and move vector
            p = cell.midpoint()
#            print("cell midpoint:", p.x())
            r0 = [p.x(),]


            # Move in POSITIVE direction into the next cell (barely)
            dr = [0.50001*self.mesh1D_dx,]
            facetExpected = 1

#            print("Calling find_facet:", "r0 =", r0, "dr =", dr, "ci =", cell.index())
            (facetCrossed, pathFraction, facetNormalVector) = dolfin_functions_solib.find_facet(self.mesh1D, r0, dr, cell.index(), returnStdArray=False)
#            print("Returned values =", facetCrossed, "pathFraction =", pathFraction, "facetNormalVector =", facetNormalVector)
            # Check facet
            self.assertEqual(facetCrossed, facetExpected, msg = "1D: facet crossed is not correct")
            # Look up the cell crossed into
            cell2Index = neighborCellIndices[facetCrossed]
#            print "ci =", cell.index(), "cell2Index =", cell2Index
            # Skip movements outside the mesh upper boundary
            if cell2Index != Mesh_C.NO_CELL:
                self.assertEqual(cell.index()+1, cell2Index, msg = "1D: cell crossed into is not correct")
                # Look up the mesh-level facet index of the facet crossed
                ifacet = self.mesh1D.cell_entity_indices_dict['facet'][cell.index()][facetCrossed]
#                print 'Cell facet', facetCrossed, 'has mesh index', ifacet

                # Check this: look up the facet (dimension = 0) shared
                # by the two cells,
                facet = df_m.MeshEntity(self.mesh1D.mesh, 0, ifacet)
                # and make sure the two cells with this facet (the
                # cells are entities having dimension 1) are the start
                # and end cells.
#                print 'facet is connected to cell', facet.entities(1)

                self.assertTrue(cell.index() in facet.entities(1) and cell2Index in facet.entities(1), msg = "1D: facet given by find_facet is not correct")
#            if cell2Index != Mesh_C.NO_CELL:ENDIF

            # Check the returned facet-normal vector against the expected vector.

            decimalPlaces = 6 # Accuracy of test
            expectedFacetNormal = [1.,]
            self.assertAlmostEqual(facetNormalVector[0], expectedFacetNormal[0], places=decimalPlaces, msg="Cell facet normal is not correct")
            numberOfFacetCrossings += 1


            # Now, move in NEGATIVE direction into the next cell (barely)
            dr = [-0.50001*self.mesh1D_dx,]
            facetExpected = 0

#            print "Calling find_facet:", "r0 =", r0, "dr =", dr, "ci =", cell.index()

            (facetCrossed, pathFraction, facetNormalVector) = dolfin_functions_solib.find_facet(self.mesh1D, r0, dr, cell.index(), returnStdArray=False)
#            print("Returned values =", facetCrossed, "pathFraction =", pathFraction, "facetNormalVector =", facetNormalVector)
            self.assertEqual(facetCrossed, facetExpected, msg = "1D facet crossed is not correct")
            # Look up the cell crossed into
            cell2Index = neighborCellIndices[facetCrossed]
            # Skip movements outside the mesh lower boundary
            if cell2Index != Mesh_C.NO_CELL:
                self.assertEqual(cell.index()-1, cell2Index, msg = "1D: cell crossed into is not correct")

                # Look up the mesh-level facet index of the facet crossed
                ifacet = self.mesh1D.cell_entity_indices_dict['facet'][cell.index()][facetCrossed]
#                print 'Cell facet', facetCrossed, 'has mesh index', ifacet

                # Check this: look up the facet (dimension = 0) shared
                # by the two cells,
                facet = df_m.MeshEntity(self.mesh1D.mesh, 0, ifacet)
                # and make sure the two cells with this facet (the
                # cells are entities having dimension 1) are the start
                # and end cells.
#                print 'facet is connected to cell', facet.entities(1)

                self.assertTrue(cell.index() in facet.entities(1) and cell2Index in facet.entities(1), msg = "1D: facet given by find_facet is not correct")

            # Check the returned facet-normal vector against the expected vector.
            decimalPlaces = 6 # Accuracy of test
            expectedFacetNormal = [-1.,]
            self.assertAlmostEqual(facetNormalVector[0], expectedFacetNormal[0], places=decimalPlaces, msg="Cell facet normal is not correct")
#            if cell2Index != Mesh_C.NO_CELL:ENDDIF
            numberOfFacetCrossings += 1
#        for cell in df_m.cells(self.mesh1D.mesh): ENDFOR

        print("\t1D mesh: tested", numberOfFacetCrossings, "facet crossings")                

        #
        # 2D: Initial point = center of cell, move = to center of neighbor cells
        #
        numberOfFacetCrossings = 0
        for cell in df_m.cells(self.mesh2D.mesh):
            neighborCellIndices = self.mesh2D.mea_object.get_cell_neighbors(cell.index())
            for nci in neighborCellIndices:
                if nci == Mesh_C.NO_CELL: # Indicates that one facet of the cell is
                                          # on a boundary, i.e., there is no
                                          # neighbor cell on that facet.
                    continue
                cellIndexExpected = nci

                # Create initial point and move vector
                p = cell.midpoint()
#                print "cell midpoint:", p.x(), p.y()
                r0 = [p.x(), p.y()]

                cell2 = self.mesh2D.cell_dict[cellIndexExpected]
                p2 = cell2.midpoint()
                dr = [p2.x()-p.x(), p2.y()-p.y()]

#                print "r0, dr=", r0, dr
#                (facetCrossed, pathFraction, facetNormalVector) = self.mesh2D.find_facet(r0, dr, cell.index())
                (facetCrossed, pathFraction, facetNormalVector) = dolfin_functions_solib.find_facet(self.mesh2D, r0, dr, cell.index())
#                print("Returned values =", facetCrossed, "pathFraction =", pathFraction, "facetNormalVector =", facetNormalVector)

#                print "start cell", cell.index(), "expected_cell", cellIndexExpected, "facet crossed", facetCrossed
                # Lookup the cell crossed from the facet number
                cell2Index = neighborCellIndices[facetCrossed]
#                print "cell2Index=", cell2Index

                self.assertEqual(cell2Index, cellIndexExpected, msg = "2D: final cell is not correct")

                # Look up the mesh-level facet index of the facet crossed
                ifacet = self.mesh2D.cell_entity_indices_dict['facet'][cell.index()][facetCrossed]
#                print 'Cell facet', facetCrossed, 'has mesh facet index', ifacet

                # Check this: look up the facet (dimension = 1) shared
                # by the two cells,
                facet = df_m.MeshEntity(self.mesh2D.mesh, 1, ifacet)
                # and make sure the two cells with this facet (the
                # cells are entities having dimension 2) are the start
                # and end cells.
#                print "Cells owning this facet:", facet.entities(2), "vs.", cell.index(), cell2Index
                self.assertTrue(cell.index() in facet.entities(2) and cell2Index in facet.entities(2), msg = "2D: facet given by find_facet is not correct")
                numberOfFacetCrossings += 1
                
        print("\t2D mesh: tested", numberOfFacetCrossings, "facet crossings")
        
        #
        # 3D: Initial point = center of cell, move = to center of neighbor cells
        #
        numberOfFacetCrossings = 0
        for cell in df_m.cells(self.mesh3D.mesh):
            neighborCellIndices = self.mesh3D.mea_object.get_cell_neighbors(cell.index())            

            # Loop on this cell's neighbors. Create a displacement
            # vector from the center of this cell to the center of the
            # neighbor.  Get the cell-index of the neighbor using (1)
            # the cell_dict[] and (2) find_facet().  These two values
            # should be the same.

            for nci in neighborCellIndices:
                if nci == Mesh_C.NO_CELL: # Indicates that one facet of the cell is
                                          # on a boundary, i.e., there is no neighbor
                                          # cell on that facet.
                    continue
                cellIndexExpected = nci

                # Create initial point and move vector
                p = cell.midpoint()
#                print "cell midpoint:", p.x(), p.y()
                r0 = [p.x(), p.y(), p.z()]

                cell2 = self.mesh3D.cell_dict[cellIndexExpected]
                p2 = cell2.midpoint()
                dr = [p2.x()-p.x(), p2.y()-p.y(), p2.z()-p.z()]

#                print "r0, dr=", r0, dr
                (facetCrossed, pathFraction, facetNormalVector) = dolfin_functions_solib.find_facet(self.mesh3D, r0, dr, cell.index())

#                print "start cell", cell.index(), "expected_cell", cellIndexExpected, "facet crossed", facetCrossed

                # Lookup the cell crossed from the facet number
                cell2Index = neighborCellIndices[facetCrossed]
#                print "cell2Index=", cell2Index

                # Check that this is the expected cell
                self.assertEqual(cell2Index, cellIndexExpected, msg = "3D: final cell is not correct")

                # Look up the mesh-level facet index of the facet crossed
#                print 'Facets of 3D mesh', self.mesh3D.cell_entity_indices_dict['facet']
#                print 'Facets of cell', self.mesh3D.cell_entity_indices_dict['facet'][cell.index()]
                ifacet = self.mesh3D.cell_entity_indices_dict['facet'][cell.index()][facetCrossed]
                # Check this: Used ifacet to look up the facet entity
                # (dimension = 2) shared by the two cells,
                facet = df_m.MeshEntity(self.mesh3D.mesh, 2, ifacet)
                # and make sure the two cells with this facet (the
                # cells are entities having dimension 3) are the start
                # and end cells.
#                print "Cells owning this facet:", facet.entities(3), "vs.", cell.index(), cell2Index
                self.assertTrue(cell.index() in facet.entities(3) and cell2Index in facet.entities(3), msg = "3D: facet given by find_facet is not correct")
                numberOfFacetCrossings += 1
                
        print("\t3D mesh: tested", numberOfFacetCrossings, "facet crossings")


        return
#    def test_2_facet_crossing(self):ENDDEF

#class TestFacetCrossing(unittest.TestCase):ENDCLASS

if __name__ == '__main__':
    unittest.main()
