#!/usr/bin/env python3

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import os
import numpy as np_m
import importlib as im_m
import unittest

import dolfin as df_m

from Dolfin_Module import Mesh_C

from UserMesh_y_Fields_FE_XYZ_Module import *

class TestFacetCrossing(unittest.TestCase):
    """Test of finding the facet crossed and desination cell due to a simulated
       displacement.

       No Particle_C object is used.
    """
    
    def setUp(self):

        # Initializations performed before each test go here...

        plotFlag = False
#        plotFlag = True
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name

        # Turn off plotting if there's no DISPLAY

        # if os.environ.get('DISPLAY') is None:
        #     plotFlag=False
        # else:
        #     plotFlag=True

        # Create 1D, 2D and 3D meshes that the particles can be tested against

        # 1d mesh input
        umi1D = UserMeshInput_C()
        pmin = -0.03
        pmax = 0.03
        cells_on_side = 4
        umi1D.pmin = df_m.Point(pmin)
        umi1D.pmax = df_m.Point(pmax)
        umi1D.cells_on_side = (cells_on_side,) # Need the comma to indicate a tuple
        self.mesh1D_dx = (pmax-pmin)/cells_on_side

        # Create mesh
        self.mesh1D = UserMesh_C(umi1D, compute_dictionaries=True, compute_cpp_arrays=False, compute_tree=True, plot_flag=plotFlag, plot_title=plotTitle + ": 1D")
#        self.mesh1D.compute_cell_vertices_dict()
#        self.mesh1D.compute_cell_dict()

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

        # Create mesh
        self.mesh2D = UserMesh_C(umi2D, compute_dictionaries=True, compute_cpp_arrays=False, compute_tree=True, plot_flag=plotFlag, plot_title=plotTitle + ": 2D")
#        self.mesh2D.compute_cell_vertices_dict()
#        self.mesh2D.compute_cell_dict()

        # 3D mesh input
        umi3D = UserMeshInput_C()
        umi3D.pmin = df_m.Point(-0.03, -0.03, -0.03)
        umi3D.pmax = df_m.Point(0.03, 0.03, 0.03)
        umi3D.cells_on_side = (4, 4, 4)
        # Create mesh
        self.mesh3D = UserMesh_C(umi3D, compute_tree=True, plot_flag=plotFlag, plot_title=plotTitle + ": 3D")
        self.mesh3D.compute_cell_dict()
        self.mesh3D.compute_cell_entity_indices_dict('vertex')
        self.mesh3D.compute_cell_entity_indices_dict('facet')
        self.mesh3D.compute_cell_facet_normals_dict() # Unit vectors normal to cell facets
        self.mesh3D.compute_cell_neighbors_dict()

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
        self.mesh1D.mesh.init(0,0)

        vertex_coords = self.mesh1D.mesh.coordinates()

        for cell in df_m.cells(self.mesh1D.mesh):
            # Obtain the facet-normal vectors for this cell
            facet_normal_vectors = self.mesh1D.cell_facet_normals_dict[cell.index()]
#            print("cell", cell.index(), "from dict, cell_facet_normals", facet_normal_vectors)

            # This returns the list of facet indices.  These are mesh-level indices, not the indices local to a cell.
#            print "cell", cell.index(), "has theses facets", self.mesh1D.cell_entity_indices_dict['facet'][cell.index()]
            mesh_dim = self.mesh1D.gdim
            for fi in range(mesh_dim+1):
                self.assertAlmostEqual(cell.normal(fi).x(), facet_normal_vectors[fi,0], msg = "1D cell facet normal x component is not correct")

#            print "cell", cell.index(), "normal0", cell.normal(0).x(), "normal1", cell.normal(1).x()

            # Print vertices of each facet
            for findx in self.mesh1D.cell_entity_indices_dict['facet'][cell.index()]:

                # A way to get the facet from it's index:
                facet = df_m.Facet(self.mesh1D.mesh, findx)
#                normal = cell.normal(findx)
                # Get the normal from the facet instead (no difference)
                normal = facet.normal()
#                print "cell", cell.index(), "facet", findx, "normal", normal.x()
                # This gives the same index:
#                print "the facet index is :",facet.index()
                # The vertices touching this facet:
                verts = facet.entities(0) # this is an ndarray, not a number of entities. DOLFIN ERROR (reported)
                
#                print "verts", verts[0]
                # Check if the facet is on a mesh boundary
#                print "exterior", facet.exterior() # Gives correct result!

#                for v in verts:
#                    print "vertex:", v, "coords:", vertex_coords[v]

#        sys.exit()

        # 2D mesh
        vertex_coords = self.mesh2D.mesh.coordinates()
        for cell in df_m.cells(self.mesh2D.mesh):
            facet_normal_vectors = self.mesh2D.cell_facet_normals_dict[cell.index()]
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
                self.assertAlmostEqual(cell.normal(fi).x(), facet_normal_vectors[fi,0], msg = "2D cell facet normal x component is not correct")
                self.assertAlmostEqual(cell.normal(fi).y(), facet_normal_vectors[fi,1], msg = "2D cell facet normal y component is not correct")

        # 3D mesh
        vertex_coords = self.mesh3D.mesh.coordinates()

        for cell in df_m.cells(self.mesh3D.mesh):
            facet_normal_vectors = self.mesh3D.cell_facet_normals_dict[cell.index()]
            verts = cell.entities(0)
            # Print values for a few cells
#            for i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,]:
            for i in [0,]:
                if cell.index() == i:
                    verts = cell.entities(0)
                    for v in verts:
#                        print "vertex:", v, "coords:", vertex_coords[v]
                        pass

                    # use the "local index" 0, 1, 2 for the triangle facets
                    # print "cell", cell.index(), "normal0", cell.normal(0).x(), cell.normal(0).y(), cell.normal(0).z()
                    # print "cell", cell.index(), "normal1", cell.normal(1).x(), cell.normal(1).y(), cell.normal(1).z()
                    # print "cell", cell.index(), "normal2", cell.normal(2).x(), cell.normal(2).y(), cell.normal(2).z()
                    # print "cell", cell.index(), "normal3", cell.normal(3).x(), cell.normal(3).y(), cell.normal(3).z()
                    pass
            # Check all the normals for all cells
            mesh_dim = self.mesh3D.gdim
            for fi in range(mesh_dim+1):
                self.assertAlmostEqual(cell.normal(fi).x(), facet_normal_vectors[fi,0], msg = "3D cell facet normal x component is not correct")
                self.assertAlmostEqual(cell.normal(fi).y(), facet_normal_vectors[fi,1], msg = "3D cell facet normal y component is not correct")
                self.assertAlmostEqual(cell.normal(fi).z(), facet_normal_vectors[fi,2], msg = "3D cell facet normal z component is not correct")

        return
#    def test_1_facet_normals(self):ENDDEF


    def test_2_facet_crossing(self):
        """ Test calculation of the cell facet crossed by a move vector.

            No Particle_C object is used.
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest:', fncName, '('+__file__+')')

        #
        # 1D: Initial point = center of cell, move = cell length
        #
        for cell in df_m.cells(self.mesh1D.mesh):
            neighbor_cell_indices = self.mesh1D.cell_neighbors_dict[cell.index()]

            # Create initial point and move vector
            p = cell.midpoint()
#            print("cell midpoint:", p.x())
            r0 = [p.x(),]

            # Move in POSITIVE direction into the next cell (barely)
            dr = [0.50001*self.mesh1D_dx,]
            facet_expected = 1

#            print("Calling find_facet:", "r0 =", r0, "dr =", dr, "ci =", cell.index())
            (facet_crossed, path_fraction, facetNormal) = self.mesh1D.find_facet(r0, dr, cell.index())
            print("Returned values =", facet_crossed, "path_fraction =", path_fraction, "facetNormal =", facetNormal)
            # Check facet
            self.assertEqual(facet_crossed, facet_expected, msg = "1D: facet crossed is not correct")
            # Look up the cell crossed into
            cell2_index = neighbor_cell_indices[facet_crossed]
#            print "ci =", cell.index(), "cell2_index =", cell2_index
            # Skip movements outside the mesh upper boundary
            if cell2_index != Mesh_C.NO_CELL:
                self.assertEqual(cell.index()+1, cell2_index, msg = "1D: cell crossed into is not correct")
                # Look up the mesh-level facet index of the facet crossed
                ifacet = self.mesh1D.cell_entity_indices_dict['facet'][cell.index()][facet_crossed]
#                print 'Cell facet', facet_crossed, 'has mesh index', ifacet

                # Check this: look up the facet (dimension = 0) shared
                # by the two cells,
                facet = df_m.MeshEntity(self.mesh1D.mesh, 0, ifacet)
                # and make sure the two cells with this facet (the
                # cells are entities having dimension 1) are the start
                # and end cells.
#                print 'facet is connected to cell', facet.entities(1)

                self.assertTrue(cell.index() in facet.entities(1) and cell2_index in facet.entities(1), msg = "1D: facet given by find_facet is not correct")
#            if cell2_index != Mesh_C.NO_CELL:ENDIF

            # Move in NEGATIVE direction into the next cell (barely)
            dr = [-0.50001*self.mesh1D_dx,]
            facet_expected = 0

#            print "Calling find_facet:", "r0 =", r0, "dr =", dr, "ci =", cell.index()
            (facet_crossed, path_fraction, facetNormal) = self.mesh1D.find_facet(r0, dr, cell.index())
            print("Returned values =", facet_crossed, "path_fraction =", path_fraction, "facetNormal =", facetNormal)            
            self.assertEqual(facet_crossed, facet_expected, msg = "1D facet crossed is not correct")
            # Look up the cell crossed into
            cell2_index = neighbor_cell_indices[facet_crossed]
            # Skip movements outside the mesh lower boundary
            if cell2_index != Mesh_C.NO_CELL:
                self.assertEqual(cell.index()-1, cell2_index, msg = "1D: cell crossed into is not correct")

                # Look up the mesh-level facet index of the facet crossed
                ifacet = self.mesh1D.cell_entity_indices_dict['facet'][cell.index()][facet_crossed]
#                print 'Cell facet', facet_crossed, 'has mesh index', ifacet

                # Check this: look up the facet (dimension = 0) shared
                # by the two cells,
                facet = df_m.MeshEntity(self.mesh1D.mesh, 0, ifacet)
                # and make sure the two cells with this facet (the
                # cells are entities having dimension 1) are the start
                # and end cells.
#                print 'facet is connected to cell', facet.entities(1)

                self.assertTrue(cell.index() in facet.entities(1) and cell2_index in facet.entities(1), msg = "1D: facet given by find_facet is not correct")
#            if cell2_index != Mesh_C.NO_CELL:ENDDIF
#        for cell in df_m.cells(self.mesh1D.mesh): ENDFOR

        #
        # 2D: Initial point = center of cell, move = to center of neighbor cells
        #
        for cell in df_m.cells(self.mesh2D.mesh):
            neighbor_cell_indices = self.mesh2D.cell_neighbors_dict[cell.index()]

            for nci in neighbor_cell_indices:
                if nci == Mesh_C.NO_CELL: # Indicates that one facet of the cell is
                                          # on a boundary, i.e., there is no
                                          # neighbor cell on that facet.
                    continue
                cell_index_expected = nci

                # Create initial point and move vector
                p = cell.midpoint()
#                print "cell midpoint:", p.x(), p.y()
                r0 = [p.x(), p.y()]

                cell2 = self.mesh2D.cell_dict[cell_index_expected]
                p2 = cell2.midpoint()
                dr = [p2.x()-p.x(), p2.y()-p.y()]

#                print "r0, dr=", r0, dr
                (facet_crossed, path_fraction, facetNormal) = self.mesh2D.find_facet(r0, dr, cell.index())

#                print "start cell", cell.index(), "expected_cell", cell_index_expected, "facet crossed", facet_crossed
                # Lookup the cell crossed from the facet number
                cell2_index = neighbor_cell_indices[facet_crossed]
#                print "cell2_index=", cell2_index

                self.assertEqual(cell2_index, cell_index_expected, msg = "2D: final cell is not correct")

                # Look up the mesh-level facet index of the facet crossed
                ifacet = self.mesh2D.cell_entity_indices_dict['facet'][cell.index()][facet_crossed]
#                print 'Cell facet', facet_crossed, 'has mesh facet index', ifacet

                # Check this: look up the facet (dimension = 1) shared
                # by the two cells,
                facet = df_m.MeshEntity(self.mesh2D.mesh, 1, ifacet)
                # and make sure the two cells with this facet (the
                # cells are entities having dimension 2) are the start
                # and end cells.
#                print "Cells owning this facet:", facet.entities(2), "vs.", cell.index(), cell2_index
                self.assertTrue(cell.index() in facet.entities(2) and cell2_index in facet.entities(2), msg = "2D: facet given by find_facet is not correct")

        #
        # 3D: Initial point = center of cell, move = to center of neighbor cells
        #
        for cell in df_m.cells(self.mesh3D.mesh):
            neighbor_cell_indices = self.mesh3D.cell_neighbors_dict[cell.index()]

            # Loop on this cell's neighbors. Create a displacement
            # vector from the center of this cell to the center of the
            # neighbor.  Get the cell-index of the neighbor using (1)
            # the cell_dict[] and (2) find_facet().  These two values
            # should be the same.

            for nci in neighbor_cell_indices:
                if nci == -1: # Indicates that one facet of the cell is
                              # on a boundary, i.e., there is no
                              # neighbor cell on that facet.
                    continue
                cell_index_expected = nci

                # Create initial point and move vector
                p = cell.midpoint()
#                print "cell midpoint:", p.x(), p.y()
                r0 = [p.x(), p.y(), p.z()]

                cell2 = self.mesh3D.cell_dict[cell_index_expected]
                p2 = cell2.midpoint()
                dr = [p2.x()-p.x(), p2.y()-p.y(), p2.z()-p.z()]

#                print "r0, dr=", r0, dr
                (facet_crossed, path_fraction, facetNormal) = self.mesh3D.find_facet(r0, dr, cell.index())

#                print "start cell", cell.index(), "expected_cell", cell_index_expected, "facet crossed", facet_crossed

                # Lookup the cell crossed from the facet number
                cell2_index = neighbor_cell_indices[facet_crossed]
#                print "cell2_index=", cell2_index

                # Check that this is the expected cell
                self.assertEqual(cell2_index, cell_index_expected, msg = "3D: final cell is not correct")

                # Look up the mesh-level facet index of the facet crossed
#                print 'Facets of 3D mesh', self.mesh3D.cell_entity_indices_dict['facet']
#                print 'Facets of cell', self.mesh3D.cell_entity_indices_dict['facet'][cell.index()]
                ifacet = self.mesh3D.cell_entity_indices_dict['facet'][cell.index()][facet_crossed]
                # Check this: Used ifacet to look up the facet entity
                # (dimension = 2) shared by the two cells,
                facet = df_m.MeshEntity(self.mesh3D.mesh, 2, ifacet)
                # and make sure the two cells with this facet (the
                # cells are entities having dimension 3) are the start
                # and end cells.
#                print "Cells owning this facet:", facet.entities(3), "vs.", cell.index(), cell2_index
                self.assertTrue(cell.index() in facet.entities(3) and cell2_index in facet.entities(3), msg = "3D: facet given by find_facet is not correct")


        return
#    def test_2_facet_crossing(self):ENDDEF

#class TestFacetCrossing(unittest.TestCase):ENDCLASS

if __name__ == '__main__':
    unittest.main()
