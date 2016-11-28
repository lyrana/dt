# Contains methods that use the DOLFIN finite-element library

#__version__ = 0.1
#__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

"""Dolfin_Module contains classes that use the DOLFIN finite-element library.
"""

import sys
import os
import dolfin as df_M
import numpy as np_M

class Mesh_C(object):
    """Mesh_C defines the mesh. UserMesh_C is a subclass.
    
       :param Mesh: a Dolfin Mesh object
       :param meshFile: The name of a file containing an XML mesh
       :param computeDictionaries: Boolean flag to compute dictionaries of mesh entities
       :param computeTree: Boolean flag to compute the mesh bounding-box tree

       :cvar gdim: Geometric dimension of the mesh
       :cvar tdim: Topological dimension of the mesh
    """

    # Static class variables

#    NO_CELL = 0xFFFFFFFF
    NO_CELL = -1
    NO_FACET = -1

    # Constructor
    def __init__(self, Mesh=None, meshFile=None, computeDictionaries=False, computeTree=False, plotFlag=False):

        # Get the mesh from the arg list
        if Mesh is not None:
            self.mesh = Mesh
        # or read the mesh from a file
        elif meshFile is not None:
            self.mesh = df_M.Mesh(meshFile)

        # This section will also execute if this __init__ is called by a
        # child class.
        if self.mesh is not None:
            if plotFlag == True:
                df_M.plot(self.mesh, axes=True)
                df_M.interactive()
            # Compute the search tree. Uses:
            #     Evaluating field probes at a point.
            #     Computing the cell index of a particle.
            if computeTree == True:
                self.bbtree = self.mesh.bounding_box_tree()

            self.gdim = self.mesh.geometry().dim()
            self.tdim = self.mesh.topology().dim()

            self.dx = np_M.empty(self.gdim, np_M.float64) # !!! This is a particle quantity requiring particle precision??  maybe using mesh precision is more consistent with the operations here since this is a mesh class. DOLFIN precision is double = float64.
            self.x0 = np_M.empty(self.gdim, np_M.float64)

            self.facet_normal_3d = np_M.empty(3, np_M.float64)

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

            # Compute lookup dictionaries
            if computeDictionaries == True:
                self.compute_cell_dict()
                self.compute_cell_entity_index_dict('vertex') # Get vertex indices from a cell index
                self.compute_cell_entity_index_dict('facet') # Get facet indices from a cell index
                self.compute_vertex_cell_dict() # Get cell indices from a vertex index

                self.compute_cell_facet_normals_dict() # Unit vectors normal to cell facets
                self.compute_cell_neighbor_dict() # Get neighbor cell indices from a cell index

        return
#    def __init__(self, Mesh=None, meshFile=None, computeDictionaries=False, computeTree=False, plotFlag=False):ENDDEF


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
    def compute_cell_dict(self):
        """Create a dictionary where the keys are cell indices, and
           the values are the corresponding cell objects.
        """

        # If the cells are numbered consecutively, this could just be an array.
        self.cell_dict = {cell.index(): cell for cell in df_M.cells(self.mesh)}

        return

#
#  Things indexed by cell number
#

    def compute_cell_entity_index_dict(self, entity_name):
        """
           Make a dictionary giving a list of the cell entity indices, indexed by
           the cell index.

           Vertex indices can be used, e.g., to obtain the x, y, z
           coordinates of the vertices of a cell.

           :param entity_name: The name of the entity: 'vertex', 'edge', 'facet'
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

        self.cell_entity_index_dict[entity_name] = dict((cell.index(), cell.entities(d)) for cell in df_M.cells(self.mesh))

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
        self.cell_vertex_dict = dict((cell.index(), cell.entities(0)) for cell in df_M.cells(self.mesh))

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
        tdim = self.tdim

        # Generate facet-cell connectivity data
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

# The following line creates a correct dictionary, but if a cell is at a boundary it
# skips that facet.  That doesn't work because we need to use the
# local facet number as an index for the list of neighbor cells.
#        self.cell_neighbor_dict = {cell.index(): sum((filter(lambda ci: ci != cell.index(), facet.entities(tdim)) for facet in df_M.facets(cell)), []) for cell in df_M.cells(self.mesh)}

        n_facets = tdim + 1
#        n_normal_coords = self.gdim

        for cell in df_M.cells(self.mesh):

            # A cell has n_facets, and each facet can have at most 1 neighbor-cell attached.
            # Use a 32-bit int for the cell indices.
            neighbor_cells = np_M.empty(n_facets, np_M.int32)

            fi = 0 # fi ranges from 0 to n_facets-1
            for facet in df_M.facets(cell):
                # Obtain the indices of cells attached to this facet
                attached_cell_indices = filter(lambda ci: ci != cell.index(), facet.entities(tdim))
                # Can only have at most 1 neighbor-cell attached
                if len(attached_cell_indices) > 0:
                    neighbor_cells[fi] = attached_cell_indices[0]
                else:
#                    if tdim ==2: print "find_facet(): tdim=", tdim, "cell_index=", cell.index(), "facet=", fi, "attached_cell_indices=", attached_cell_indices
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
        self.vertex_cell_dict = dict((v.index(), [c.index() for c in df_M.cells(v)]) for v in df_M.vertices(self.mesh))

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
        
        for cell in df_M.cells(self.mesh):
            normals = np_M.empty(shape=(n_facets, n_normal_coords), dtype = np_M.float64)
            for fi in range(n_facets):
                facet_normal_3d[0] = cell.normal(fi).x()
                facet_normal_3d[1] = cell.normal(fi).y()
                facet_normal_3d[2] = cell.normal(fi).z()
#                normals[fi][:] = facet_normal_3d
                normals[fi] = facet_normal_3d[0:n_normal_coords]
            self.cell_facet_normals_dict[cell.index()] = normals

#        self.cell_facet_normals = dict((cell.index(), [cell.normal(fi) for fi in range(self.tdim+1)]) for cell in df_M.cells(self.mesh))

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
            p = df_M.Point(point['x'], point['y'], point['z'])
        elif self.gdim == 2:
            p = df_M.Point(point['x'], point['y'])
        else:
            p = df_M.Point(point['x'])

#        first = self.bbtree.compute_first_entity_collision(df_M.Point(p)) # not compute_first_collision()
        first = self.bbtree.compute_first_entity_collision(p) # not compute_first_collision()

        fncname = sys._getframe().f_code.co_name + '():'
#        print fncname, "cell index =", first

        # The maximum value of a 32-bit int is 0xFFFFFFFF (-1)
        # (eight Fs and F = 1111)
#        if first == 0xFFFFFFFF:
        if first == Mesh_C.NO_CELL:
            print "fncname: particle is out of bounds"
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

        fncname = sys._getframe().f_code.co_name + '():'

        for ip in xrange(points.shape[0]):
        # pseg[i] is a 'record' containing the 'x', 'y', 'z', 'vx', 'vy',... values of ith ithem
        # so pseg[i][0:3] is 'x', 'y', 'z'.

            p = [points[ip][d] for d in range(gdim)]
            ploc = Point(p)
            compute_cell_index(ploc)

        first = self.bbtree.compute_first_entity_collision(point)
        print fncname, "cell index =", first

        return
        
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
            p = df_M.Point(point['x'], point['y'], point['z'])
        elif self.gdim == 2:
            p = df_M.Point(point['x'], point['y'])
        else:
            p = df_M.Point(point['x'])

#        return cell.contains(df_M.Point(self.gdim, point))
        return cell.contains(df_M.Point(p))

#class Mesh_C(object):
    def find_facet(self, r0, dr, cell_index):
        """Compute the cell facet crossed in traveling along a
           displacement vector dr from position r0.

           :param r0: initial position, as a tuple of coordinates (x),
                      (x,y), or (x, y, z) with field names 'x', 'y',
                      'z'

           :param dr: move vector, as a tuple of coordinates (x),
                      (x,y), or (x, y, z) with field names 'x', 'y',
                      'z'
           
           :param cell_index: the (global?) index of the cell that the
                              particle started in.

        :returns: number of facet traversed, or None,
                  fraction of the path that's in this cell.
        :rtype: int, float
        """

        fncname = sys._getframe().f_code.co_name + '():'

        dim = self.gdim

        dx = self.dx
        x0 = self.x0
        dx[:] = dr[0:dim] # Convert displacement vector to a vector with
                   # dimension equal to the number of mesh dimensions.
        x0[:] = r0[0:dim]

        facet = Mesh_C.NO_FACET

        path_fraction = 1.0
        distance_to_facet = 0.0

        if np_M.dot(dx,dx) == 0.0:
            return facet, path_fraction

        cell = self.cell_dict[cell_index]

        # The coordinates of this cell's vertices, in (x, y, z) tuple form
        vertex_coords = cell.get_vertex_coordinates().reshape((-1, dim))

#        print "find_facet(): vertex_coords=", vertex_coords

        # The facet-normals of this cell:
        facet_normal_vectors = self.cell_facet_normals_dict[cell_index]        

        if dim == 3:

            # 3D mesh: There are 4 facets, and the normal vector is 3D.

            # Test if the plane of facet 0 is crossed:
            n0_dot_dx = np_M.dot(facet_normal_vectors[0], dx)
#            print "find_facet(): facet_normal_vectors[0] = ", facet_normal_vectors[0], "dx=", dx, "n0_dot_dx=", n0_dot_dx
            if n0_dot_dx > 0:
                # The vector from the starting point to a vertex in the facet plane
                vec_to_facet = np_M.subtract(vertex_coords[1], x0)
                # The normal distance to the facet plane
                distance_to_facet = np_M.dot(facet_normal_vectors[0], vec_to_facet)
#                print "f0 find_facet(): vec_to_facet=", vec_to_facet, "distance_to_facet=", distance_to_facet
                if distance_to_facet < path_fraction*n0_dot_dx:
                    facet = 0
                    path_fraction = distance_to_facet/n0_dot_dx
#                    print "f0 find_facet(): facet=", facet, "path_fraction=", path_fraction

            # Next, test if the plane of facet 1 is crossed:
            n1_dot_dx = np_M.dot(facet_normal_vectors[1], dx)
#            print "find_facet(): n1_dot_dx=", n1_dot_dx
            if n1_dot_dx > 0:
                # The vector from the starting point to a vertex in the facet plane
                vec_to_facet = np_M.subtract(vertex_coords[2], x0)
                # The normal distance to the facet plane
                distance_to_facet = np_M.dot(facet_normal_vectors[1], vec_to_facet)
                if distance_to_facet < path_fraction*n1_dot_dx:
                    facet = 1
                    path_fraction = distance_to_facet/n1_dot_dx

            # Next, test if the plane of facet 2 is crossed:
            n2_dot_dx = np_M.dot(facet_normal_vectors[2], dx)
#            print "find_facet(): n2_dot_dx=", n2_dot_dx
            if n2_dot_dx > 0:
                # The vector from the starting point to a vertex in the facet plane
                vec_to_facet = np_M.subtract(vertex_coords[3], x0)
                # The normal distance to the facet plane
                distance_to_facet = np_M.dot(facet_normal_vectors[2], vec_to_facet)
                if distance_to_facet < path_fraction*n2_dot_dx:
                    facet = 2
                    path_fraction = distance_to_facet/n2_dot_dx

            # Next, test if the plane of facet 3 is crossed:
            n3_dot_dx = np_M.dot(facet_normal_vectors[3], dx)
#            print "find_facet(): n3_dot_dx=", n3_dot_dx
            if n3_dot_dx > 0:
                # The vector from the starting point to a vertex in the facet plane
                vec_to_facet = np_M.subtract(vertex_coords[0], x0)
                # The normal distance to the facet plane
                distance_to_facet = np_M.dot(facet_normal_vectors[3], vec_to_facet)
                if distance_to_facet < path_fraction*n3_dot_dx:
                    facet = 3
                    path_fraction = distance_to_facet/n3_dot_dx

            if path_fraction > 1.0:
                facet = Mesh_C.NO_FACET
                print fncname, "!!! path_fraction > 1.0:", path_fraction
                sys.exit()
            return facet, path_fraction

        elif dim == 2:

            # 2D mesh: There are 3 facets, and the normal vector is 2D.

            # Test if the plane of facet 0 is crossed:
            n0_dot_dx = np_M.dot(facet_normal_vectors[0], dx)
#            print "find_facet(): facet_normal_vectors[0] = ", facet_normal_vectors[0], "dx=", dx, "n0_dot_dx=", n0_dot_dx
            if n0_dot_dx > 0:
                # The vector from the starting point to a vertex in the facet plane
                vec_to_facet = np_M.subtract(vertex_coords[1], x0)
                # The normal distance to the facet plane
                distance_to_facet = np_M.dot(facet_normal_vectors[0], vec_to_facet)
#                print "f0 find_facet(): vec_to_facet=", vec_to_facet, "distance_to_facet=", distance_to_facet
                if distance_to_facet < path_fraction*n0_dot_dx:
                    facet = 0
                    path_fraction = distance_to_facet/n0_dot_dx
#                    print "f0 find_facet(): facet=", facet, "path_fraction=", path_fraction

            # Next, test if the plane of facet 1 is crossed:
            n1_dot_dx = np_M.dot(facet_normal_vectors[1], dx)
#            print "find_facet(): n1_dot_dx=", n1_dot_dx
            if n1_dot_dx > 0:
                # The vector from the starting point to a vertex in the facet plane
                vec_to_facet = np_M.subtract(vertex_coords[2], x0)
                # The normal distance to the facet plane
                distance_to_facet = np_M.dot(facet_normal_vectors[1], vec_to_facet)
                if distance_to_facet < path_fraction*n1_dot_dx:
                    facet = 1
                    path_fraction = distance_to_facet/n1_dot_dx

            # Next, test if the plane of facet 2 is crossed:
            n2_dot_dx = np_M.dot(facet_normal_vectors[2], dx)
#            print "find_facet(): n2_dot_dx=", n2_dot_dx
            if n2_dot_dx > 0:
                # The vector from the starting point to a vertex in the facet plane
                vec_to_facet = np_M.subtract(vertex_coords[0], x0)
                # The normal distance to the facet plane
                distance_to_facet = np_M.dot(facet_normal_vectors[2], vec_to_facet)
                if distance_to_facet < path_fraction*n2_dot_dx:
                    facet = 2
                    path_fraction = distance_to_facet/n2_dot_dx

            if path_fraction > 1.0:
                facet = Mesh_C.NO_FACET
                print fncname, "!!! path_fraction > 1.0:", path_fraction
                sys.exit()
            return facet, path_fraction
        else:

            # 1D mesh: There are 2 facets, and the normal vector is 1D.

            # Test if facet 0 is crossed:
            n0_dot_dx = facet_normal_vectors[0]*dx[0]
#            print "find_facet()", "n0_dot_dx=", n0_dot_dx
            if n0_dot_dx > 0:
                # The vector from the starting point to a vertex in the facet plane
                vec_to_facet = vertex_coords[0] - x0[0]
                # The normal distance to the facet plane
                distance_to_facet = facet_normal_vectors[0]*vec_to_facet
                if distance_to_facet < path_fraction*n0_dot_dx:
                    facet=0
                    path_fraction = distance_to_facet/n0_dot_dx
                    return facet, path_fraction

            # Next, test if facet 1 is crossed:
            n1_dot_dx = facet_normal_vectors[1]*dx[0]
#            print "find_facet()", "n1_dot_dx=", n1_dot_dx
            if n1_dot_dx > 0:
                # The vector from the starting point to a vertex in the facet plane
                vec_to_facet = vertex_coords[1] - x0[0]
                # The normal distance to the facet plane
                distance_to_facet = facet_normal_vectors[1]*vec_to_facet
#                print "find_facet()", "distance_to_facet:", distance_to_facet, "vec_to_facet:", vec_to_facet, "facet_normal_vectors[1]:", facet_normal_vectors[1]
                if distance_to_facet < path_fraction*n1_dot_dx:
                    facet=1
                    path_fraction = distance_to_facet/n1_dot_dx
#                    print "find_facet", "n1_dot_dx=", n1_dot_dx, "vec_to_facet=", vec_to_facet, "distance_to_facet=", distance_to_facet
                    return facet, path_fraction

            if path_fraction > 1.0:
                facet = Mesh_C.NO_FACET
                print fncname, "!!! path_fraction > 1.0:", path_fraction
                sys.exit()
            return facet, path_fraction
#    def find_facet(self, r0, dr, cell_index):ENDDEF

#class Mesh_C(object):END


class Field_C(object):
    """Field_C uses the DOLFIN library to represent
       fields defined on finite-element spaces.

       :param meshCI: A Mesh_C object.  Contains the finite-element
                      mesh where the field is defined.
       :param element_type: The name of the basis function family.
       :param element_degree: The order of the polynomial in the basis
                              functions.
       :param field_type: 'scalar' or 'vector'.

    """

    def __init__(self, meshCI, element_type, element_degree, field_type):

        self.meshCI = meshCI
        self.mesh_gdim = meshCI.mesh.geometry().dim()
        self.mesh_tdim = meshCI.mesh.topology().dim()

        if field_type is 'scalar':
            self.function_space = df_M.FunctionSpace(meshCI.mesh, element_type, element_degree)
        elif field_type is 'vector':
            self.function_space = df_M.VectorFunctionSpace(meshCI.mesh, element_type, element_degree)

        # Create the finite-element function for this field
        self.function = df_M.Function(self.function_space)
        self.function_values = self.function.vector()

# Q: Is it a good idea to set the following variables?
# A: Then you don't have to know Dolfin function name when writing DT scripts.
        self.element_type = element_type
        self.element_degree = element_degree

# dofmap here?        

        return
#    def __init__(self, meshCI, element_type, element_degree, field_type):END

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
        fieldValue = np_M.empty(len(field_at_points.dtype.fields), dtype=field_at_points.dtype[0])
#        fieldValue = np_M.empty(self.mesh_gdim, dtype=vFpoints.dtype[0])

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
        """This function takes the inner product of a delta-function
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
            point = df_M.Point(p['x'], p['y'], p['z'])
        elif self.mesh_gdim == 2:
            point = df_M.Point(p['x'], p['y'])
        else:
            point = df_M.Point(p['x'])

        df_M.PointSource(self.function_space, point, p['weight']).apply(self.function_values)

        # get cell from the particle index
#        c = self.mesh3DCI.cell_dict[pseg[ip]['cell_index']]
        # Evaluate the basis at the location

        return
#    def integrate_delta_function(self, p):END

#class Field_C(object):END


class PoissonSolve_C(object):
    """PoissonSolve_C uses the DOLFIN library for solving field
       equations with the finite-element method.  The methods here are
       ones that the user doesn't have to modify for set up different
       simulations.  (Functions that the user can modify are in the
       UserPoissonSolve_C class in the UserMesh module.)
       
       Methods:
           solve_for_phi()
           compute_negE
    """

    def __init__(self):

        return

#class PoissonSolve_C(object):
    def solve_for_phi(self, plotFlag=False):
        """Solve Poisson's equation for the electric potential.
        """

# Note that "scalarFunctionSpace" must be defined in the child class
#        self.phi = df_M.Function(self.scalarFunctionSpace)
#        df_M.solve(self.a == self.L, self.phi, self.bcs)

# Set the RHS (charge-density sources)
        if self.charge_density is None:
            f = df_M.Constant(0.0)
        else:
            f = self.charge_density
        self.L = f*self.v*df_M.dx

# Compute solution
        df_M.solve(self.a == self.L, self.u, self.bcs, solver_parameters=self.solver_parameters)

# Plot the electric potential
        if plotFlag == True:
            df_M.plot(self.u)
            df_M.interactive()

# Compute -E
        # This is supposed to test if neg_electric_field has been allocated.
        if self.neg_electric_field is not None:
            self.compute_negE(plotFlag)

# Plot radial component of -E:
#        negE = df_M.grad(u)

# Project onto a piece-wise constant-in-element function
#        negE2 = df_M.project(E, df_M.VectorFunctionSpace(mesh, 'DG', 0))

#        negE2_x = E2[0]
#        negE2_y = E2[1]

#        df_M.plot(negE2_x, title="negE2_x")
#        df_M.plot(negE2_y, title="negE2_y")

#        df_M.interactive()
        return

#class PoissonSolve_C(object):
    def compute_negE(self, plotFlag=False):
        """negE is the gradient of the electrostatic potential.
        """

        negE = df_M.grad(self.u)

        function_space = self.neg_electric_field.function_space
        self.neg_electric_field.function = df_M.project(negE, function_space)
        if plotFlag == True:
            df_M.plot(self.neg_electric_field.function)
            df_M.interactive()

        return self.neg_electric_field.function
# how do you change sign?: see Epoints

#class PoissonSolve_C(object):END
