# Contains methods that use the DOLFIN finite-element library

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
# Alphabetical list of exports
__all__ = ['Field_C',
           'Mesh_C',
           'PoissonSolve_C',
           'RectangularRegion_C',
           'SphericalRegion_C',
           'WholeMesh_C']

"""Dolfin_Module contains classes that use the DOLFIN finite-element library.
"""

import sys
import os
import dolfin as df_m
import numpy as np_m

#STARTCLASS
class Mesh_C(object):
    """Mesh_C defines the mesh. UserMesh_C is a subclass that adds
       boundary-conditions, etc.
    
    """

    # Static class variables

#    NO_CELL = 0xFFFFFFFF
    NO_CELL = -1
    NO_FACET = -1

    # Constructor

# Examine the need for the Mesh and mesh_to_copy args:
    def __init__(self, Mesh=None, mesh_to_copy=None, mesh_file=None, field_boundary_marker_file=None, compute_dictionaries=False, compute_tree=False, plot_flag=False, plot_title="Mesh"):
        """
           :param Mesh: a Dolfin Mesh object
           :param mesh_file: The name of a file containing an XML mesh
           :param compute_dictionaries: Boolean flag to compute dictionaries of mesh entities
           :param compute_tree: Boolean flag to compute the mesh bounding-box tree

           :cvar gdim: Geometric dimension of the mesh
           :cvar tdim: Topological dimension of the mesh
        """

        # Get a ref to the mesh from the arglist,
        if Mesh is not None:
            self.mesh = Mesh
        # or make a copy of the arglist mesh,
        elif mesh_to_copy is not None:
            self.mesh = df_m.Mesh(mesh_to_copy)
        # or read the mesh from a file.
        elif mesh_file is not None:
            self.mesh = df_m.Mesh(mesh_file)

        # This section will also execute if this __init__ is called by a
        # child class.
        if self.mesh is not None:
            if plot_flag is True:
                df_m.plot(self.mesh, title=plot_title, axes=True)
                df_m.interactive()
            # Compute the search tree. Uses:
            #     Evaluating field probes at a point.
            #     Computing the cell index of a particle.
            if compute_tree is True:
                self.bbtree = self.mesh.bounding_box_tree()

            self.gdim = self.mesh.geometry().dim()
            self.tdim = self.mesh.topology().dim()

            self.dx = np_m.empty(self.gdim, np_m.float64) # !!! This is a particle quantity requiring particle precision??  maybe using mesh precision is more consistent with the operations here since this is a mesh class. DOLFIN precision is double = float64.
            self.x0 = np_m.empty(self.gdim, np_m.float64)

            self.facet_normal_3d = np_m.empty(3, np_m.float64)

            # Don't need to declare all of these class variables here
            # provided they're created by the dict() function in the
            # class methods.  They're listed here only for
            # documentation.
            self.cell_dict = {}

            self.entity_dimension = {'vertex': 0, 'edge': 1, 'facet': self.tdim-1}
            self.cell_entity_index_dict = {}
            self.vertex_cell_dict = {}
            self.cell_facet_normals_dict = {}
            self.cell_neighbor_dict = {}

            self.cell_volume_dict = {}

            # Compute lookup dictionaries
            if compute_dictionaries is True:
                self.compute_cell_dict()
                self.compute_cell_entity_index_dict('vertex') # Get vertex indices from a cell index
                self.compute_cell_entity_index_dict('facet') # Get facet indices from a cell index
                self.compute_vertex_cell_dict() # Get cell indices from a vertex index

                self.compute_cell_facet_normals_dict() # Unit vectors normal to cell facets
                self.compute_cell_neighbor_dict() # Get neighbor cell indices from a cell index

                self.compute_cell_volume_dict() # Compute cell volumes indexed by cell index

            # Read the field boundary marker from a file, if one is provided
            if field_boundary_marker_file is not None:
                fieldBoundaryMarkerFile_df = df_m.File(field_boundary_marker_file)
                fieldBoundaryMarker = df_m.MeshFunction('size_t', self.mesh, self.mesh.topology().dim()-1)
                fieldBoundaryMarkerFile_df >> fieldBoundaryMarker
                self.field_boundary_marker = fieldBoundaryMarker

        return
#    def __init__(self, Mesh=None, mesh_file=None, compute_dictionaries=False, compute_tree=False, plot_flag=False):ENDDEF


 # Compute the search tree for this mesh

        
#   BoundingBoxTree tree;
#             self.pmesh_tree = self.pmesh.bounding_box_tree()

#   tree.build(mesh, dim); // leave out the dim argument for default == cells

# You can then use the tree to compute different kinds of collisions:

#   compute_collisions(const Point& point)
#   compute_entity_collisions(const Point& point, const Mesh& mesh)
#   compute_first_collision(const Point& point) const
#   compute_first_entity_collision(const Point& point, const Mesh& mesh)
#   compute_closest_entity(const Point& point, const Mesh& mesh)
#   compute_closest_point(const Point& point)

# In your case, I assume you need to compute at least one (the first)
# collision between a point and the cells of the mesh:

#   unsigned int cell_index = tree.compute_first_entity_collision(point, mesh);

#class Mesh_C(object):
    def __str__(self):
        """Print the class members.
        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():'
        print fncName

        print_string = ' mesh_tdim = ' + str(self.mesh_tdim)
        print_string += '\nmesh_gdim = ' + str(self.mesh_gdim)

        return print_string
#     def __str__(self):ENDDEF

    def copy(self):
        """ The 'copy constructor'.
        """

        return copy.deepcopy(self)

#class Mesh_C(object):
    def compute_cell_dict(self):
        """Create a dictionary where the keys are cell mesh-indices, and
           the values are the corresponding cell objects.
        """

        # If the cells are numbered consecutively, this could just be an array.
        self.cell_dict = {cell.index(): cell for cell in df_m.cells(self.mesh)}

        return
#    def compute_cell_dict(self):ENDDEF

#
#  Things indexed by cell number
#

    def compute_cell_entity_index_dict(self, entity_name):
        """
           Make a dictionary where the keys are cell indices and the values are
           lists of the cell entity indices.

           Vertex indices can be used, e.g., to obtain the x, y, z
           coordinates of the vertices of a cell.

           :param entity_name: The name of the cell entity that will be in the
                               value list, e.g., 'vertex', 'edge', 'facet'
           :type entity_name: string

        """

        tdim = self.tdim
        d = self.entity_dimension[entity_name]

        # Initialize the necessary connectivity
        self.mesh.init(d, tdim)

        # For each cell, get a list of the entities in it
        # (cell.entities(d) returns a list of indices).
        # A vertex has topological dimension 0.
        # A facet has topological codimension d-1.

        # Type issue: Had to convert numpy.uint32 to long (using tolist()) and
        # then to int (using map()).
        # This is because a FacetFunction wants the facet index as an int.
        self.cell_entity_index_dict[entity_name] = dict((cell.index(), map(int, cell.entities(d).tolist())) for cell in df_m.cells(self.mesh))

        # E.g., use the dictionary to print the vertex coordinates:

        # Can't this be done with builtin dolfin iterators?

        # coordinates = self.mesh.coordinates()
        # for c in range(self.mesh.num_cells()):
        #     print "Vertices of cell", c
        #     for vertex in self.cell_vertex_dict[c]:
        #         print "\t", vertex, coordinates[vertex]

        return
#    def compute_cell_entity_index_dict(self, entity_name):ENDDEF

#class Mesh_C(object):
    def compute_cell_vertex_dict(self):
        """
           Make a dictionary giving a list of the cell vertex indices, indexed by
           the cell index.  The vertex indices can be used, e.g., to
           obtain the x, y, z coordinates of the vertices.
        """

        # For each cell, get a list of the vertices in it
        # A vertex has topological dimension 0
        self.cell_vertex_dict = dict((cell.index(), cell.entities(0)) for cell in df_m.cells(self.mesh))

        # Use the dictionary to print the vertex coordinates:

        # coordinates = self.mesh.coordinates()
        # for c in range(self.mesh.num_cells()):
        #     print "Vertices of cell", c
        #     for vertex in self.cell_vertex_dict[c]:
        #         print "\t", vertex, coordinates[vertex]

        return
#    def compute_cell_vertex_dict(self):ENDDEF

#class Mesh_C(object):
    def compute_cell_neighbor_dict(self):
        """Make a dictionary that gives the indices of cells that
           share a common facet.

           The dictionary is indexed by cell index.
        
        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        tdim = self.tdim

        # Generate facet-to-cell connectivity data
        self.mesh.init(tdim - 1, tdim)

        # For every cell, build a list of cells that are connected to its facets
        # but are not the iterated cell.

        # Create a dictionary, where the keys are the cell indices, and the values are a list of the cell indices of touching cells.
        # Double loop: outer one over cells in the mesh, inner one over facets of the cell.
        # cells(mesh) are the cells in the mesh. 
        # cell.index() is the local index of the cell.
        # facets(cell) are the facets of a cell.
        # facet.entities(tdim) are the entities of dimension tdim touching that facet

        # Breakdown:
        #   cell.index() is the key in the dictionary.

        #   (lambda ci: ci != cell.index(): 
        #       A function that tests if ci is not the current cell.

        #   filter(test, list to test): Apply the above test to the
        #       indices of dimension tdim that are attached to
        #       'facet'.  Return a list of the items that pass the
        #       test.

        #   sum(list, []): Concatenate the lists obtained for each
        #       cell facet, removing extra levels of brackets. These
        #       are the values in the dictionary. E.g., [1,2], [3,4,5]
        #       becomes [1, 2, 3, 4, 5].
        #  


# Used "tdim" for quantities independent of geometry, like indices.

# The following commented-out line creates a correct dictionary, but if a
# cell is at a boundary it skips that facet.  That doesn't work for
# our purpose because we need to use the local facet number as an
# index for the list of neighbor cells.

#        self.cell_neighbor_dict = {cell.index(): sum((filter(lambda ci: ci != cell.index(), facet.entities(tdim)) for facet in df_m.facets(cell)), []) for cell in df_m.cells(self.mesh)}

        n_facets = tdim + 1
#        n_normal_coords = self.gdim

        for cell in df_m.cells(self.mesh):

            # Make a numpy array to hold the indices of the neighbor cells.
            # A cell has n_facets, and each facet can have at most 1 neighbor-cell attached.
            # Use a 32-bit int for the cell indices.
            neighbor_cells = np_m.empty(n_facets, np_m.int32)

            fi = 0 # fi ranges from 0 to n_facets-1
            for facet in df_m.facets(cell):
                # Obtain the indices of cells attached to this facet
                # See breakdown of this line above.
                attached_cell_indices = filter(lambda ci: ci != cell.index(), facet.entities(tdim))
                # Can only have at most 1 neighbor-cell attached
                if len(attached_cell_indices) > 0:
                    neighbor_cells[fi] = attached_cell_indices[0]
                else:
#                    if tdim ==2: print fncName, tdim, "cell_index=", cell.index(), "facet=", fi, "attached_cell_indices=", attached_cell_indices
#                    neighbor_cells[fi] = -1
                    neighbor_cells[fi] = Mesh_C.NO_CELL
                fi += 1

            self.cell_neighbor_dict[cell.index()] = neighbor_cells

#        print cell_neighbor_dict

        return
#    def compute_cell_neighbor_dict(self):ENDDEF

#
#  Things indexed by vertex number
#

#class Mesh_C(object):
    def compute_vertex_cell_dict(self):
        """
           Make a dictionary giving a list of the indices of cells
           that have a given vertex index.
        """

        # For each vertex, get a list of the cells that contain this vertex
        self.vertex_cell_dict = dict((v.index(), [c.index() for c in df_m.cells(v)]) for v in df_m.vertices(self.mesh))

        return

#class Mesh_C(object):
    def compute_cell_facet_normals_dict(self):
        """
           Make a dictionary giving a list of cell facet-normal
           vectors, indexed by the cell index.

           NO:The normals are unit vectors represented by Points.  

           The normals are unit vectors with the number of coordinates
           equal to the number of coordinates of the mesh.

           The dictionary looks like: {cell_index: [n0, n1, ...]}
           where n0, n1, ... are the unit normal vectors on each facet.

        """

        # For each cell, get a list of the vertices in it
        # (cell.entities(d) returns a list of indices).
        # A vertex has dimension 0.

        # cell.normal(fi) returns a Point representing the unit normal to the fi-th facet of the cell

        n_facets = self.gdim+1
        n_normal_coords = self.gdim

        facet_normal_3d = self.facet_normal_3d
        
        for cell in df_m.cells(self.mesh):
            normals = np_m.empty(shape=(n_facets, n_normal_coords), dtype = np_m.float64)
            for fi in range(n_facets):
                facet_normal_3d[0] = cell.normal(fi).x()
                facet_normal_3d[1] = cell.normal(fi).y()
                facet_normal_3d[2] = cell.normal(fi).z()
#                normals[fi][:] = facet_normal_3d
                normals[fi] = facet_normal_3d[0:n_normal_coords]
            self.cell_facet_normals_dict[cell.index()] = normals

#        self.cell_facet_normals = dict((cell.index(), [cell.normal(fi) for fi in range(self.tdim+1)]) for cell in df_m.cells(self.mesh))

        return
#    def compute_cell_facet_normals_dict(self):ENDDEF

#
#  Functions of arbitrary point positions
#
#class Mesh_C(object):
    def compute_cell_index(self, point):
        """Given a coordinate tuple with field names 'x', 'y', 'z',
           find the (local) index of the cell that contains it.  If
           the point lies on the boundary between two cells, the index
           of the first one found is returned.

           :param point: (x), (x,y), or (x,y,z) coordinate tuple

           :returns: cell index containing the point, or NO_CELL if no cell
                     on the mesh contains it.

        """

        if self.gdim == 3:
            p = df_m.Point(point['x'], point['y'], point['z'])
        elif self.gdim == 2:
            p = df_m.Point(point['x'], point['y'])
        else:
            p = df_m.Point(point['x'])

#        first = self.bbtree.compute_first_entity_collision(df_m.Point(p)) # not compute_first_collision()
        first = self.bbtree.compute_first_entity_collision(p) # not compute_first_collision()

        fncName = sys._getframe().f_code.co_name + '():'
#        print fncName, "cell index =", first

        # The maximum value of a 32-bit int is 0xFFFFFFFF (-1)
        # (eight Fs and F = 1111)
#        if first == 0xFFFFFFFF:
        if first == Mesh_C.NO_CELL:
            print fncName, ": particle is out of bounds"
            sys.exit()

        return first
#    def compute_cell_index(self, point):ENDDEF

#class Mesh_C(object):
    def compute_cell_indices(self, points):
        """Given an array of points, find the (local) index of the cell that
           contains each one.  If the point lies on the boundary between two
           cells, the index of the first cell found is returned.
           Args: point is of dolfin type Point.
        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        for ip in xrange(points.shape[0]):
        # pseg[i] is a 'record' containing the 'x', 'y', 'z', 'vx', 'vy',... values of ith ithem
        # so pseg[i][0:3] is 'x', 'y', 'z'.

            # Doesn't work!
            p = [points[ip][d] for d in range(gdim)]
            ploc = Point(p)
            compute_cell_index(ploc)

        first = self.bbtree.compute_first_entity_collision(point)
        print fncName, "cell index =", first

        return

#class Mesh_C(object):
    def compute_cell_volume_dict(self):
        """
           Make a dictionary giving the cell volume, indexed by
           the cell index.
        """

        self.cell_volume_dict = dict((cell.index(), cell.volume()) for cell in df_m.cells(self.mesh))

        return
#    def compute_cell_volume_dict(self):ENDDEF

    def is_inside(self, point, cell_index):
        """Check if a point lies within a give cell

        :param point: a tuple of coordinates (x), (x,y), or (x, y, z)
                      with field names 'x', 'y', 'z'
        :param cell_index: a (global?) cell index

        :returns: True iff the point is inside the cell.
        :rtype: bool
        """

        cell = self.cell_dict[cell_index]

        if self.gdim == 3:
            p = df_m.Point(point['x'], point['y'], point['z'])
        elif self.gdim == 2:
            p = df_m.Point(point['x'], point['y'])
        else:
            p = df_m.Point(point['x'])

#        return cell.contains(df_m.Point(self.gdim, point))
        return cell.contains(df_m.Point(p))

#class Mesh_C(object):
    def find_facet(self, r0, dr, cell_index):
        """Compute the cell facet crossed in traveling along a
           displacement vector dr from position r0.

           The facet crossed is the one with the smallest value of the
           fraction (normal distance to facet plane)/(total distance
           traveled normal to facet plane)

           1. Compute the distance a particle moves in the directions
           of a facet (Lf).

           2. If it moved closer to it, compute the distance to the
           facet (Df).

           3. If Df < Lf, the particle crossed the plane that the
           facet is in.

           4. If the particle crosses more than one facet plane, the
           facet crossed is the one with the smallest ratio Df/Lf,
           since it reached that facet first.

           :param r0: initial position, as a tuple of coordinates (x),
                      (x,y), or (x, y, z)

           :param dr: move vector, as a tuple of coordinates (x),
                      (x,y), or (x, y, z)
           
           :param cell_index: the (global?) index of the cell that the
                              particle started in.

        :returns: number of facet traversed, or None,
                  fraction of the path that's in this cell.
        :rtype: int, float
        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        gDim = self.gdim # The geometric dimension of the mesh

        dx = self.dx
        x0 = self.x0
        dx[:] = dr[0:gDim] # Convert displacement vector to a vector with
                          # dimension equal to the number of mesh dimensions.
        x0[:] = r0[0:gDim]
#        x0[:] = r0[gDim:2*gDim]

        facet = Mesh_C.NO_FACET

        pathFraction = 1.0
        distanceToFacet = 0.0

        if np_m.dot(dx,dx) == 0.0:
            return facet, pathFraction

        cell = self.cell_dict[cell_index]

        # The coordinates of all the vertices in this cell, in (x, y,
        # z) tuple form
        vertex_coords = cell.get_vertex_coordinates().reshape((-1, gDim))

#        print "find_facet(): vertex_coords=", vertex_coords

        # The facet-normals of this cell:
        facet_normal_vectors = self.cell_facet_normals_dict[cell_index]        

        if gDim == 3:

            # 3D mesh: There are 4 facets indexed 0,1,2,3, and the normal vector is 3D.

            # Test if the plane of facet 0 is crossed:
            n0_dot_dx = np_m.dot(facet_normal_vectors[0], dx)
#            print "find_facet(): facet_normal_vectors[0] = ", facet_normal_vectors[0], "dx=", dx, "n0_dot_dx=", n0_dot_dx
            if n0_dot_dx > 0:
                # The vector from the starting point to a vertex in the facet plane
                vecToFacet = np_m.subtract(vertex_coords[1], x0)
                # The normal distance to the facet plane
                distanceToFacet = np_m.dot(facet_normal_vectors[0], vecToFacet)
#                print "f0 find_facet(): vecToFacet=", vecToFacet, "distanceToFacet=", distanceToFacet
                if distanceToFacet < 0.0:
                    print fncName, "!!! Bad value for distanceToFacet:", distanceToFacet, "flipping the sign!!!"
                    # Assume this is round-off error and flip the sign
                    distanceToFacet = -distanceToFacet
                if distanceToFacet < pathFraction*n0_dot_dx:
                    facet = 0
                    pathFraction = distanceToFacet/n0_dot_dx
#                    print "f0 find_facet(): facet=", facet, "pathFraction=", pathFraction

            # Next, test if the plane of facet 1 is crossed:
            n1_dot_dx = np_m.dot(facet_normal_vectors[1], dx)
#            print "find_facet(): n1_dot_dx=", n1_dot_dx
            if n1_dot_dx > 0:
                # The vector from the starting point to a vertex in the facet plane
                vecToFacet = np_m.subtract(vertex_coords[2], x0)
                # The normal distance to the facet plane
                distanceToFacet = np_m.dot(facet_normal_vectors[1], vecToFacet)
                if distanceToFacet < 0.0:
                    print fncName, "!!! Bad value for distanceToFacet:", distanceToFacet, "flipping the sign!!!"
                    # Assume this is round-off error and flip the sign
                    distanceToFacet = -distanceToFacet
                if distanceToFacet < pathFraction*n1_dot_dx:
                    facet = 1
                    pathFraction = distanceToFacet/n1_dot_dx

            # Next, test if the plane of facet 2 is crossed:
            n2_dot_dx = np_m.dot(facet_normal_vectors[2], dx)
#            print "find_facet(): n2_dot_dx=", n2_dot_dx
            if n2_dot_dx > 0:
                # The vector from the starting point to a vertex in the facet plane
                vecToFacet = np_m.subtract(vertex_coords[3], x0)
                # The normal distance to the facet plane
                distanceToFacet = np_m.dot(facet_normal_vectors[2], vecToFacet)
                if distanceToFacet < 0.0:
                    print fncName, "!!! Bad value for distanceToFacet:", distanceToFacet, "flipping the sign!!!"
                    # Assume this is round-off error and flip the sign
                    distanceToFacet = -distanceToFacet
                if distanceToFacet < pathFraction*n2_dot_dx:
                    facet = 2
                    pathFraction = distanceToFacet/n2_dot_dx

            # Next, test if the plane of facet 3 is crossed:
            n3_dot_dx = np_m.dot(facet_normal_vectors[3], dx)
#            print "find_facet(): n3_dot_dx=", n3_dot_dx
            if n3_dot_dx > 0:
                # The vector from the starting point to a vertex in the facet plane
                vecToFacet = np_m.subtract(vertex_coords[0], x0)
                # The normal distance to the facet plane
                distanceToFacet = np_m.dot(facet_normal_vectors[3], vecToFacet)
                if distanceToFacet < 0.0:
                    print fncName, "!!! Bad value for distanceToFacet:", distanceToFacet, "flipping the sign!!!"
                    # Assume this is round-off error and flip the sign
                    distanceToFacet = -distanceToFacet
                if distanceToFacet < pathFraction*n3_dot_dx:
                    facet = 3
                    pathFraction = distanceToFacet/n3_dot_dx

            if pathFraction > 1.0 or pathFraction < 0.0:
                print fncName, "!!! Bad value for pathFraction:", pathFraction
                sys.exit()
            return facet, pathFraction

        elif gDim == 2:

            # 2D mesh: There are 3 facets indexed 0,1,2, (opposite the
            # vertices with these same indices) and the normal vector is
            # 2D.

            # Test if the plane of facet 0 is crossed.
            # Distance traveled normal to facet 0:
            n0_dot_dx = np_m.dot(facet_normal_vectors[0], dx)
#            print "find_facet(): facet_normal_vectors[0] = ", facet_normal_vectors[0], "dx=", dx, "n0_dot_dx=", n0_dot_dx
            if n0_dot_dx > 0:
                # The vector from the starting point to a vertex in the facet plane
                vecToFacet = np_m.subtract(vertex_coords[1], x0)
                # The normal distance to the facet plane
                distanceToFacet = np_m.dot(facet_normal_vectors[0], vecToFacet)
                if distanceToFacet < 0.0:
                    print fncName, "!!! Bad value for distanceToFacet:", distanceToFacet, "flipping the sign!!!"
                    # Assume this is round-off error and flip the sign
                    distanceToFacet = -distanceToFacet
#                print "f0 find_facet(): vecToFacet=", vecToFacet, "distanceToFacet=", distanceToFacet
                if distanceToFacet < pathFraction*n0_dot_dx:
                    facet = 0
                    # Compute fraction of the path to this facet in this cell:
                    pathFraction = distanceToFacet/n0_dot_dx
#                    print "f0 find_facet(): facet=", facet, "pathFraction=", pathFraction

            # Next, test if the plane of facet 1 is crossed.
            # Distance traveled normal to facet 0:
            n1_dot_dx = np_m.dot(facet_normal_vectors[1], dx)
#            print "find_facet(): n1_dot_dx=", n1_dot_dx
            if n1_dot_dx > 0:
                # The vector from the starting point to a vertex in the facet plane
                vecToFacet = np_m.subtract(vertex_coords[2], x0)
                # The normal distance to the facet plane
                distanceToFacet = np_m.dot(facet_normal_vectors[1], vecToFacet)
                # If the path-fraction to this facet's plane is
                # smaller than for facet 0, this may be the one crossed.
                if distanceToFacet < 0.0:
                    print fncName, "!!! Bad value for distanceToFacet:", distanceToFacet, "flipping the sign!!!"
                    # Assume this is round-off error and flip the sign
                    distanceToFacet = -distanceToFacet
                if distanceToFacet < pathFraction*n1_dot_dx:
                    facet = 1
                    # Compute fraction of the path to this facet in this cell:
                    pathFraction = distanceToFacet/n1_dot_dx

            # Next, test if the plane of facet 2 is crossed.
            # Distance traveled normal to facet 2:
            n2_dot_dx = np_m.dot(facet_normal_vectors[2], dx)
#            print "find_facet(): n2_dot_dx=", n2_dot_dx
            if n2_dot_dx > 0:
                # The vector from the starting point to a vertex in the facet plane
                vecToFacet = np_m.subtract(vertex_coords[0], x0)
                # The normal distance to the facet plane
                distanceToFacet = np_m.dot(facet_normal_vectors[2], vecToFacet)
                if distanceToFacet < 0.0:
                    print fncName, "!!! Bad value for distanceToFacet:", distanceToFacet, "flipping the sign!!!"
                    # Assume this is round-off error and flip the sign
                    distanceToFacet = -distanceToFacet
                # If the path-fraction to this facet's plane is
                # smaller than for facet 0 and 1, this is the one crossed.
                if distanceToFacet < pathFraction*n2_dot_dx:
                    facet = 2
                    pathFraction = distanceToFacet/n2_dot_dx

            if pathFraction > 1.0 or pathFraction < 0.0:
                print fncName, "!!! Bad value for pathFraction:", pathFraction, "!!!"
                sys.exit()
#                if abs(pathFraction) < TINY_PATH_FRACTION:
#                    pathFraction = 0.0 # this will trigger a check on the cell index?
#                else:
#                    sys.exit()
            return facet, pathFraction
        else:

            # 1D mesh: There are 2 facets indexed 0 1, and the normal vector is 1D.

            # Test if facet 0 is crossed:
            n0_dot_dx = facet_normal_vectors[0]*dx[0]
#            print "find_facet()", "n0_dot_dx=", n0_dot_dx
            if n0_dot_dx > 0:
                # The vector from the starting point to a vertex in the facet plane
                vecToFacet = vertex_coords[0] - x0[0]
                # The normal distance to the facet plane
                distanceToFacet = facet_normal_vectors[0]*vecToFacet
                if distanceToFacet < 0.0:
                    print fncName, "!!! Bad value for distanceToFacet:", distanceToFacet, "flipping the sign!!!"
                    # Assume this is round-off error and flip the sign
                    distanceToFacet = -distanceToFacet
                if distanceToFacet < pathFraction*n0_dot_dx:
                    facet=0
                    pathFraction = distanceToFacet/n0_dot_dx
                    return facet, pathFraction

            # Next, test if facet 1 is crossed:
            n1_dot_dx = facet_normal_vectors[1]*dx[0]
#            print "find_facet()", "n1_dot_dx=", n1_dot_dx
            if n1_dot_dx > 0:
                # The vector from the starting point to a vertex in the facet plane
                vecToFacet = vertex_coords[1] - x0[0]
                # The normal distance to the facet plane
                distanceToFacet = facet_normal_vectors[1]*vecToFacet
                if distanceToFacet < 0.0:
                    print fncName, "!!! Bad value for distanceToFacet:", distanceToFacet, "flipping the sign!!!"
                    # Assume this is round-off error and flip the sign
                    distanceToFacet = -distanceToFacet
#                print "find_facet()", "distanceToFacet:", distanceToFacet, "vecToFacet:", vecToFacet, "facet_normal_vectors[1]:", facet_normal_vectors[1]
                if distanceToFacet < pathFraction*n1_dot_dx:
                    facet=1
                    pathFraction = distanceToFacet/n1_dot_dx
#                    print "find_facet", "n1_dot_dx=", n1_dot_dx, "vecToFacet=", vecToFacet, "distanceToFacet=", distanceToFacet
                    return facet, pathFraction

            if pathFraction > 1.0 or pathFraction < 0.0:
                print fncName, "!!! Bad value for pathFraction:", pathFraction
                sys.exit()
            return facet, pathFraction
#    def find_facet(self, r0, dr, cell_index):ENDDEF

#class Mesh_C(object):ENDCLASS


#STARTCLASS
class Field_C(object):
    """Field_C uses the DOLFIN library to represent
       fields defined on finite-element spaces.

       :param mesh_M: A Mesh_C object.  Contains the finite-element
                      mesh where the field is defined.
       :param element_type: The name of the basis function family.
       :param element_degree: The order of the polynomial in the basis
                              functions.
       :param field_type: 'scalar' or 'vector'.

    """

    def __init__(self, mesh_M, element_type, element_degree, field_type):

        self.mesh_M = mesh_M
        self.mesh_gdim = mesh_M.mesh.geometry().dim()
        self.mesh_tdim = mesh_M.mesh.topology().dim()

        if field_type is 'scalar':
            self.function_space = df_m.FunctionSpace(self.mesh_M.mesh, element_type, element_degree)
        elif field_type is 'vector':
            self.function_space = df_m.VectorFunctionSpace(self.mesh_M.mesh, element_type, element_degree)
        self.field_type = field_type

        # Create the finite-element function for this field
        self.function = df_m.Function(self.function_space)
        self.function_values = self.function.vector()

# Q: Is it a good idea to set the following variables?
# A: Then you don't have to know Dolfin function name when writing DT scripts.
        self.element_type = element_type
        self.element_degree = element_degree

# dofmap here?        

        return
#    def __init__(self, mesh_M, element_type, element_degree, field_type):END

#class Field_C(object):
    def __str__(self):
        """Print the class members.
        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():'
        print fncName

        print_string = 'mesh_tdim = ' + str(self.mesh_tdim)
        print_string += '\nmesh_gdim = ' + str(self.mesh_gdim)
        print_string += '\nfield_type = ' + str(self.field_type)
        print_string += '\nelement_type = ' + str(self.element_type)
        print_string += '\nelement_degree = ' + str(self.element_degree)

        return print_string
#     def __str__(self):ENDDEF

#class Field_C(object):
    def set_to_value(self, value):
        """Set all the field values to the given value

           :param value: Value to set the field to.
        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        # Assignment to an FE function: FBook p. 25:
        self.function.vector()[:] = value

        return
#    def set_to_value(self, value):ENDDEF

#class Field_C(object):
    def multiply_add(self, field_to_add, multiplier):
        """Add multiplier*field_to_add to the field values in the current object.

            E.g., this is used to obtain charge-density on a mesh from number-density
            on the mesh.

           :param field_to_add: A Field_C object to be added to the current Field_C
                                object.
           :param multiplier: Scalar value that multiplies the field_to_add values
        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        # Assignment to a Dolfin mesh function: FBook p. 25:
        self.function.vector()[:] += multiplier*field_to_add.function.vector()

        return
#    def multiply_add(self, value):ENDDEF

#    def interpolate_vectorField_to_points(self, points, vFpoints):
#class Field_C(object):
    def interpolate_field_to_points(self, points, field_at_points):
        """Interpolate fields to given points.

        :param points: Array of points.
        :type points: ndarray
        :param field_at_points: The calculated field values at the points.
        :type field_at_points: ndarray

        :returns: None
           
        """

        # Temporary ndarray for the interpolated field at a point
        # Length is number of components of field (1 for scalar; 2 for 2D vector, etc.) 
        fieldValue = np_m.empty(len(field_at_points.dtype.fields), dtype=field_at_points.dtype[0])
#        fieldValue = np_m.empty(self.mesh_gdim, dtype=vFpoints.dtype[0])

# Loop on the points:

        for ip in xrange(points.shape[0]):

            # Get the coordinates to locate the point on the mesh.
            p = [points[ip][d] for d in range(self.mesh_gdim)]

            # Evaluate the finite-element function at the point p
            self.function(p, values=fieldValue)
            field_at_points[ip] = fieldValue
#        print 'field_at_points = ', field_at_points
#        return field_at_points[0:npoints]
#        return field_at_points
        return
#    def interpolate_field_to_points(self, points, field_at_points):END

#class Field_C(object):
    def integrate_delta_function(self, p):
        """This function takes the 'inner product' of a delta-function
           located at point p and the basis functions used by the
           field object.  This occurs, e.g., in computing the source
           density for Poisson's equations, where the basis functions
           are the "test" functions. The field must be a scalar.

           :param p: A point inside the domain of the function.
           :type p: An object with data fields 'x', ('y', 'z'),
                    'weight', e.g., a particle.

           :returns: None.  The DOFs of the field in the cell
                     containing the point are incremented by the value of the
                     basis function evaluated at the point, multiplied by the
                     point's weight.

        """

        # Should be able to use p['cell_index'] here instead of having
        # to search for the cell that contains the point.
        if self.mesh_gdim == 3:
            point = df_m.Point(p['x'], p['y'], p['z'])
        elif self.mesh_gdim == 2:
            point = df_m.Point(p['x'], p['y'])
        else:
            point = df_m.Point(p['x'])

        df_m.PointSource(self.function_space, point, p['weight']).apply(self.function_values)

        # get cell from the particle index
#        c = self.mesh3DCI.cell_dict[pseg[ip]['cell_index']]
        # Evaluate the basis at the location

        return
#    def integrate_delta_function(self, p):END


#class Field_C(object):
    def interpolate_delta_function_to_dofs(self, p):
        """This function takes the 'inner product' of a delta-function located
           at point p and the basis functions used by the field object.

           This can be used e.g., to assemble the space-charge source vector
           {rho*u_i*dx} for the variational form of Poisson's equation from a
           distribution of delta-function particles.  The basis functions {u_i}
           are the "test" functions. The field rho is a scalar.

           For non-Cartesian coordinates, the integral to obtain the inner product
           will contain a non-unity Jacobian factor in order to obtain the correct
           volume element.

           :param p: A point inside the domain of the function.
           :type p: An object with data fields 'x', ('y', 'z'),
                    'weight', e.g., a particle.

           :cvar Jhat: XX A factor that comes from integrating the delta function times
                       the volume element. The value is non-unity for non-Cartesian
                       coordinates.  Factors like 2\pi and 4\pi are dropped if they
                       appear on both sides of the field PDE. They should be put back
                       to get the correct physical values of density, etc.

           :returns: None.  The DOFs of the field in the cell
                     containing the point are incremented by the value of the
                     basis function evaluated at the point, multiplied by the
                     point's weight.

        """

# Another version of this function would loop on cells instead of looping on particles
# Then coords wouldn't have to be retrieved multiple times

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        cellIndex = p['cell_index']
        el = self.function_space.element()
        cell = df_m.Cell(self.mesh_M.mesh, cellIndex)
# QQQ: why is get_vertex_coordinates() used instead of get_coordinate_dofs()?
# AAA: The DoFs are the same as the vertices for CG1 elements
#      I tried both and they give the same output for CG1 elements.
#      
#        coordinate_dofs = cell.get_vertex_coordinates()
        coordinate_dofs = cell.get_coordinate_dofs()

        # Compute the basis-function values at the given point, p.
        # "basis_values" stores the values of the basis_functions at p.
        # el.space_dimension() is the number of basis functions

        nBF = el.space_dimension()
        # A numpy array is needed to pass this value to evaluate_basis_all() below.
        numberOfBasisFunctions = np_m.zeros(1, dtype=np_m.intc)
        numberOfBasisFunctions[0] = nBF

        # Make scratch for basisValues
        basisValues = np_m.zeros(nBF, dtype=float)
        # Make scratch for point
        point = np_m.zeros(self.mesh_gdim, dtype=float)
        for i in range(self.mesh_gdim):
            point[i] = p[i]

        w = p['weight']

        # Compute the "normalized" Jacobian factor
        # if self.mesh_M.coordinate_system == 'Cartesian':
        #     Jhat = 1.0
        # elif self.mesh_M.coordinate_system == '1D-spherical-radius':
        #     radius = p['x']
        #     Jhat = radius*radius # Drop the 4\pi factor.  Put it back where physical
        #                       # density values are needed.

        el.evaluate_basis_all(basisValues, point, coordinate_dofs, cell.orientation()) 
#        print fncName, basisValues

        # Now insert the values for this cell into the global density array.

        # Follow C++ source in PointSource.cpp:

        # Get the dof indices, i.e., the indices where the dof values of
        # densities from this cell are located.
        # Note: in general, these are not the numbered the same as the vertices

        dofIndices = self.function_space.dofmap().cell_dofs(cellIndex) # type: numpy.ndarray
        # This doesnt work, because array() contains a *copy* of vector(), not
        # the underlying data.
#        for i in range(len(dofIndices)):
#            self.function.vector().array()[dofIndices[i]] += basisValues[i]
# This works:
        # Sum the densities into these array locations
#        self.function.vector()[dofIndices] += w*basisValues*Jhat
        self.function.vector()[dofIndices] += w*basisValues

# Is something additional needed in parallel? .apply("add")?

# Cant get this to work; complains of a type mismatch:
#        self.function.vector().add_local(basisValues, numberOfBasisFunctions, dofIndices)

        return
#    def interpolate_delta_function_to_dofs(self, p):ENDDEF


#class Field_C(object):ENDCLASS


#STARTCLASS
class PoissonSolve_C(object):
    """PoissonSolve_C uses the DOLFIN library for solving field equations with
       the finite-element method.  The methods in this parent class are
       ones that the user doesn't have to modify for set up different
       simulations.  (Functions that the user can modify are in the
       UserPoissonSolve_C class in the UserMesh_... module.)
       
       Methods:
           solve_for_phi()
           compute_negE

       Notes on the finite-element solution method:

       The variational form of the field equation is:
                  a(u,v) = L(v)
       a(u,v) is a bilinear form in the functions 'u' and 'v'. It contains
       known coefficients that may be functions of space and time.
       'u' is the unknown solution (the "trial function").
       L(v) is a linear form in 'v', where 'v' is any function in the "test
       function" space.

       The above variational form is turned into a matrix equation
                  Ax = b
       by using finite-dimensional trial and test-function spaces [1].
       The trial-function space is spanned by the basis functions {u[i]}.
       The test-function space is spanned by the basis functions {v[j]}.

       The matrix equation Ax = b denotes a system of equations where the
       i'th equation is
                  A[i,j]*x[j] = b[i] (summation on j)

       Note that the term A[i,j] is a(u[j], v[i]). (see [2], for order and [3]
       for definition of a(u,v)), so that the i'th row of A, i.e., A[i,*], is
       obtained from a(u, v[i]).

       'x' is the column-vector of unknown coefficients used to express 'u' as a
       linear sum of the basis functions of the trial space, {u[j]}.  Thus,
       x[j] is the coefficient of basis function u[j], i,e., u = x[j]*u[j],
       summed over j.

       b[i] is obtained by evaluating L(v[i]) (which involves an integral over
       the domain [4]).

       Making the matrix A[i,j] from 'a', {u[j]}, and {v[i]}, and making the
       vector b[i] from 'L' and {v[i]} is called "assembly".

       References:
       [1] FBOOK, p. 309 (Sec. 17.4)
       [2] FBOOK, p. 140 (Sec. 5.7)
       [3] FBOOK, p. 4 (Eq. 1.10)
       [4] FBOOK, p. 4 (Eq. 1.11)

    """

    def __init__(self):
        """The problem-specific initialization is done in the child classes.

        """

        # If the bilinear form a(u,v) has time-independent coefficients, the
        # matrix A only needs to be assembled once, which redundant work.

        # A[i,j] is formed from {u[j]}, {v[i]}, and any coefficients in
        # a(u,v). Note the order of i, j here.
        if self.pde_has_constant_coeffs is True:
            self.A = df_m.assemble(self.a)
            for bc in self.bcs:
                bc.apply(self.A)

        return
#    def __init__(self, pde_has_constant_coeffs=False):ENDDEF

#class PoissonSolve_C(object):
    def set_constant_source(self, charge_density_value):
        """Set the source term in the field equation to a constant value.

           A source has to be of type 'Expression', which is a Dolfin class [1].

           References:
           [1] FBOOK, p. 196 (Sec. 103.6).

        """

        source_expression = df_m.Constant(charge_density)
        return source_expression*self.v*df_m.dx
#    def set_constant_source(self):ENDDEF

#class PoissonSolve_C(object):
#    def set_source_function(self, charge_density_expression):
    def assemble_source_expression(self, source_density):
        """Create a source vector for the RHS of the field equation, given a
           constant density, or a string expression giving the density as a
           function of the spatial coordinates.

           The density is first turned into an 'Expression' (a Dolfin class
           [1]).  Then it is used to create a linear form, namely, the integral
           of the density times a function belonging to the test-function space.
           This linear form is applied to the basis functions of the
           test-function space to obtain the source vector for the variational
           form of the field equation.  The latter process is called 'assembly'
           of the source vector.

           An example of an Expression is [2]:
           f = df_m.Expression("sin(x[0])*cos(x[1]")

           Time-dependent Expressions can also be defined [2].

           An example of a linear form is:
           f*v*dx
           where v is in the test-function space.

           References:
           [1] FBOOK, p. 196 (Sec. 10.3.6).
           [2] FBOOK, p. 197 (Sec. 10.3.6).

           :param source_density: A constant or a string expressing a function of
                                  the spatial coordinates {x[0], x[1], x[2]}

           :returns: The linear form corresponding to the charge density.

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        # Check the type of 'source_density'
        if source_density is None:
            source_expression = df_m.Constant(0.0)
        elif isinstance(source_density, str):
            source_expression = df_m.Expression(source_density)
        elif isinstance(source_density, float):
            source_expression = df_m.Constant(source_density)
        elif isinstance(source_density, int):
            source_expression = df_m.Constant(float(source_density))
        else:
            errorMsg = "%s Parameter source_density is of unexpected type" % fncName
            sys.exit(errorMsg)
        
        L = source_expression*self.v*df_m.dx
        self.b = df_m.assemble(L)

        if source_density is not None:
            for bc in self.bcs:
                bc.apply(self.b)

#        return b
        return
#    def assemble_source_expression(self, source_density):ENDDEF

#class PoissonSolve_C(object):
    def solve_for_phi(self, plot_flag=False, plot_title="solve_for_phi"):
        """Solve Poisson's equation for the electric potential.

           The variational form of the equation is:
                      a(u,v) = L(v)
           a(u,v) is a bilinear form in 'u' and 'v'
           'u' is the unknown solution ("trial function")
           L(v) is a linear form in 'v'
           v is any function in the "test function" space.

           The variational form is turned into a matrix equation:
                      Ax = b
           by using finite-dimensional trial and test-function spaces.

           Making the matrix 'A' from 'a' and the vector 'b' from 'L' is called
           "assembly".

           :cvar U: An alias for the solution vector.

        """

        ### Assemble the A matrix for the LHS ###

        # If a(u,v) has coefficients that change, then it needs to be assembled
        # again. If the coefficients don't change, then the one-time assembly
        # takes place when UserPoissonSolve... is initialized.
        if self.pde_has_constant_coeffs is False:
            self.A = df_m.assemble(self.a)
#            print "A is of type", type(self.A)
            for bc in self.bcs:
                bc.apply(self.A)

        ### Assemble the b vector for the RHS ###

        # Assembled charge sources (e.g., from delta-function kinetic particles)
        if self.assembled_charge is not None:
            # Create a Vector and initialize it with Function.vector()
            self.b = df_m.Vector(self.assembled_charge.vector())
            # Apply the boundary conditions
            for bc in self.bcs:
                bc.apply(self.b)

        # U is the solution vector
        U = self.u.vector()

        ### Solve the Laplace/Poisson equation ###

        df_m.solve(self.A, U, self.b)

# Plot the electric potential
        if plot_flag is True:
            df_m.plot(self.u, title=plot_title+": phi")
            df_m.interactive()

# Compute -E
        # This is supposed to test if neg_electric_field has been allocated.
        if self.neg_electric_field is not None:
            self.compute_negE(plot_flag, plot_title)

# Plot radial component of -E:
#        negE = df_m.grad(u)

# Project onto a piece-wise constant-in-element function
#        negE2 = df_m.project(E, df_m.VectorFunctionSpace(mesh, 'DG', 0))

#        negE2_x = E2[0]
#        negE2_y = E2[1]

#        df_m.plot(negE2_x, title="negE2_x")
#        df_m.plot(negE2_y, title="negE2_y")

#        df_m.interactive()
        return

#class PoissonSolve_C(object):
    def compute_negE(self, plot_flag=False, plot_title="compute_negE"):
        """negE is the gradient of the electrostatic potential.

           :cvar negE: The negative gradient of the electrostatic potential

        """

        negE = df_m.grad(self.u)

        function_space = self.neg_electric_field.function_space
        self.neg_electric_field.function = df_m.project(negE, function_space)
        if plot_flag is True:
            df_m.plot(self.neg_electric_field.function, title=plot_title+": -E")
            df_m.interactive()

        return self.neg_electric_field.function
# how do you change sign?: see Epoints

#class PoissonSolve_C(object):ENDCLASS

#STARTCLASS
class FacetSubDomain_C(df_m.SubDomain):
    """
       :param Mesh:

       :cvar gdim: The number of coordinates needed to specify a point in the
                   underlying configuration space

    """

    # Static class variables

    # The SUBDOMAIN_INDEX is incremented for each new FacetSubDomain_C
    FACET_SUBDOMAIN_INDEX = 0

#class FacetSubDomain_C(df_m.SubDomain):ENDCLASS

#STARTCLASS
class CellSet_C(object):
    """A CellSet_C produces and stores a set of mesh cells.

       A cell set is defined using a function called "inside()", which
       is provided by the subclasses of CellSet_C.  The "inside()" function
       is used to determine if a given cell is included in the set.  If
       "inside()" is not defined in the subclass, the cell set is the *whole
       mesh*.

       Examples of CellSet_C subclasses are regions bounded by planes, or
       spherical regions.
    """

    # Static class variables

    # The SUBDOMAIN_INDEX is incremented for each new CellSet_C
#needed?    CELL_SUBDOMAIN_INDEX = 0

    def __init__(self, mesh_M, plot_flag=False):
        """Create a CellSet_C from the given mesh.
           :param mesh_M: A Mesh_C object.
           :param plot_flag: Boolean flag to plot the subdomain

           :cvar cell_list: An iterable object containing a set of cells.
           :cvar gdim: The number of coordinates needed to specify a point in the
                       underlying configuration space
           :cvar ncell: Number of cells in the list

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():\n'

        plotFlag = plot_flag

        # Use gdim, e.g., for the x,y,z coordinates of a point
        self.gdim = mesh_M.gdim

        # If an inside() function is not defined, the cell list is the entire mesh.
        if hasattr(self, 'inside'):
            self.cell_list = []
            for cell in df_m.cells(mesh_M.mesh):
                p = cell.midpoint()
                if self.inside(p):      
                    self.cell_list.append(cell)
        else:
            self.cell_list = [cell for cell in df_m.cells(mesh_M.mesh)] # The whole mesh

        self.ncell = len(self.cell_list)
        # Check that there's at least 1 cell in the list
        if self.ncell == 0:
            errorMsg = "(DnT ERROR) There are no cells in self.cell_list"
            raise RuntimeError(errorMsg)

# mid-point, circumradius, etc.
        # Allocate storage
        self.cell_index = np_m.empty(self.ncell, dtype = np_m.int32)
        self.midpoint = np_m.empty(shape=(self.ncell, self.gdim), dtype = np_m.float64)
        self.volume = np_m.empty(self.ncell, dtype = np_m.float64)
        self.radius = np_m.empty(self.ncell, dtype = np_m.float64)

        # Note that icell is just a counter, not the Dolfin cell index.
        # One could instead create an Interator over the cells in the set.
        for icell in xrange(self.ncell):
            c = self.cell_list[icell]
            c_index = c.index()
            self.cell_index[icell] = c_index
# 12aug17 note: now have precomputed vol for each cell in mesh
#            self.volume[icell] = c.volume()
            self.volume[icell] = mesh_M.cell_volume_dict[c_index]
            midp = c.midpoint()
            for i in range(self.gdim):
                self.midpoint[icell][i] = midp[i]
            self.radius[icell] = c.circumradius()
            
        return
#    def __init__(self, mesh, plot_flag=False):ENDDEF

    def in_cell(self, point, index):
        """Test if point is in cell[index].

           :param point: A numpy array containing x, y, z, or a subset
           :param index: The index into self.cell_list list

           :cvar cell: A Dolfin Cell object

           :returns: True iff point is in cell cell_index
           :rtype: bool

        """

        cell = self.cell_list[index]

        return cell.contains(df_m.Point(point))
#    def in_cell(self, point, cell_index):ENDDEF

#class CellSet_C(object):ENDCLASS


### Useful standard shapes for defining mesh regions.

#STARTCLASS
class RectangularRegion_C(CellSet_C):
    """The RectangularRegion_C class is a specialized CellSet_C for specifying
       (e.g.) cells that generate particles.

       This works in 1D, 2D, 3D

    """

    def __init__(self, mesh_M, pmin, pmax, print_flag=False):
        """Initialize a source with point or planar boundaries.

           :param float pmin: The lower corner of the source region. This is a scalar for 1D, and a tuple for 2D, 3D
           :param float pmax: The upper corner of the source region. This is a scalar for 1D, and a tuple for 2D, 3D

        """

        # Use gdim, e.g., for the x,y,z coordinates of a point
        self.gdim = mesh_M.gdim

        self.pmin = pmin
        self.pmax = pmax

        plotFlag = print_flag

        # Call the base class constructor
        super(self.__class__, self).__init__(mesh_M, plot_flag=plotFlag)

        return
#    def __init__(self, mesh_M, xmin, xmax):ENDDEF

#class RectangularRegion_C(CellSet_C):
    def inside(self, x):
        """Used to define what cells are in this subdomain.

           :param x: The point to be checked

           :returns: True if x[] is inside the subdomain

       """
        tol = 1.0e-10

        if self.gdim == 3:
            return x[0] >= self.pmin[0] - tol and x[0] <= self.pmax[0] + tol \
               and x[1] >= self.pmin[1] - tol and x[1] <= self.pmax[1] + tol \
               and x[2] >= self.pmin[2] - tol and x[2] <= self.pmax[2] + tol
        elif self.gdim == 2:
            return x[0] >= self.pmin[0] - tol and x[0] <= self.pmax[0] + tol \
               and x[1] >= self.pmin[1] - tol and x[1] <= self.pmax[1] + tol
        else:
            return x[0] >= self.pmin - tol and x[0] <= self.pmax + tol

        # The source is in the interval [xmin, xmax]
# def inside(self, x):ENDDEF

#class RectangularRegion_C(CellSet_C):ENDCLASS
    def __str__(self):
        """Print the class members.

           :cvar gdim: Number of coordinates specifying a point in Cartesian space
           :cvar ncell: Number of cells in the list

        """

        fncName = '('+__file__+') ' + self.__class__.__name__ + "." + sys._getframe().f_code.co_name + '():' # No \n to end this since it's printed on its own line
        print fncName

        print_string = ' gdim = ' + str(self.gdim)
        print_string += '\n pmin = ' + str(self.pmin)
        print_string += '\n pmax = ' + str(self.pmax)
        print_string += '\n ncell = ' + str(self.ncell)
        print_string += '\n cell_list(size) = ' + str(type(self.cell_list)) + '('+str(len(self.cell_list))+')'
        print_string += '\n volume(size) = ' + str(type(self.volume)) + '('+str(len(self.cell_list))+')'
        print_string += '\n midpoint(size) = ' + str(type(self.midpoint)) + '('+str(len(self.cell_list))+')'

        return print_string
#     def __str__(self):ENDDEF

#class RectangularRegion_C(CellSet_C):ENDCLASS

#STARTCLASS
class SphericalRegion_C(CellSet_C):
    """
       Define a source with a circular or spherical boundary.

       This should be generalized to 3D

    """

    def __init__(self, mesh_M, center, radius, plot_flag=False):

        self.center = center
        self.radius = radius
        plotFlag = plot_flag

        # Call the base class constructor
        super(self.__class__, self).__init__(mesh_M, plot_flag=plotFlag)

        return
#    def __init__(self, center, radius):ENDDEF

#class SphericalRegion_C(CellSet_C):
    def inside(self, x):
        """Used to define what cells are in this region.

           :param x: The point to be checked

           :returns: True if x[] is inside the region

       """
        tol = 1.0e-10

        # The source is in the interval [xmin, xmax]
        return x[0] >= self.center - self.radius - tol and x[0] <= self.center + self.radius + tol
#    def inside(self, x):ENDDEF

#class SphericalRegion_C(CellSet_C):ENDCLASS

#STARTCLASS
class WholeMesh_C(CellSet_C):
    """
       Define a cell region that's the whole mesh
    """

    def __init__(self, mesh_M):

        # Call the base class constructor
        super(self.__class__, self).__init__(mesh_M)

        return
#    def __init__(self, center, radius):ENDDEF

#class WholeMesh_C(CellSet_C):ENDCLASS


#STARTCLASS
class CellSubDomain_C2(df_m.SubDomain):
    """A CellSubDomain_C is a subdomain of the mesh made up of cells.
       
       A Dolfin SubDomain is used to define a spatial region by implementing 

           (1) an inside() function to test points x[] as being in or not in the
           region.

           (2) a mark() function that puts values (size_t, int, double, bool)
           onto entities belonging to the SubDomain.  The values are applied
           using a Dolfin MeshFunction defined on the mesh entities.  The
           MeshFunction stores the marker values.

           An example is to use integer markers on cells to identify cell
           subdomains.  These markers can be used in making functions over
           finite-element spaces that depend on the subdomains (FBook,
           p. 62). E.g., thermal conductivity that varies from region to region.

       Here, we use a Dolfin CellFunction to mark cells, and then make a list of
       the marked cells.

    """

    # Static class variables

    # The SUBDOMAIN_INDEX is incremented for each new CellSubDomain_C
    CELL_SUBDOMAIN_INDEX = 0

    def __init__(self, mesh, plot_flag=False):
        """Create a CellSubDomain_C from the given mesh.
           :param mesh: A Mesh_C object.
           :param plot_flag: Boolean flag to plot the subdomain

        """

        for cell in df_m.cells(mesh):
            p = cell.midpoint()

        ## Create a CellFunction to mark the SubDomain
        cellFunction = df_m.CellFunction('size_t', mesh)

        # Set default value to zero
        cellFunction.set_all(0)

        # Mark the cells in the SubDomain with a 1
        self.mark(cellFunction, 1)

        # Create a list of the cells in this SubDomain

        for icell in xrange(len(cellFunction.array())):
            if cellFunction.array()[icell] == 1:
                pass

#        self.cell_function = cellFunction

        return
#    def __init__(self, mesh, plot_flag=False):ENDDEF

#class CellSubDomain_C2(df_m.SubDomain):ENDCLASS
