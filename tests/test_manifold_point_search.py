#!/usr/bin/env python
import numpy
from dolfin import *

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

def test_manifold_point_search():
    # Setup simple two-cell surface in 3d where the cells have different orientations
    cellname = "triangle"
    tdim = 2
    gdim = 3
    vertices = [
        (0.0, 0.0, 1.0),
        (1.0, 1.0, 1.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        ]
    cells = [
        (0, 1, 2),
        (0, 1, 3),
        ]

    mesh = Mesh()
    me = MeshEditor()
    me.open(mesh, cellname, tdim, gdim)
    me.init_vertices(len(vertices))
    for i, v in enumerate(vertices):
        me.add_vertex(i, *v)
    me.init_cells(len(cells))
    for i, c in enumerate(cells):
        me.add_cell(i, *c)
    me.close()
    mesh.init_cell_orientations(Expression(("0.0", "0.0", "1.0")))

#    plot(mesh, title='test mesh', axes=True)
#    interactive()

    midpoints = [
        (2.0/3.0, 1.0/3.0, 2.0/3.0), # midpoint of cell 0, works ok
        (1.0/3.0, 2.0/3.0, 2.0/3.0), # midpoint of cell 1, getting v(p) = ( 0.55555556  0.44444444  0.88888889 ) here
        ]
    bb = mesh.bounding_box_tree()
    for expected_cellid in (0, 1):
        point = Point(numpy.asarray(midpoints[expected_cellid]))
        found_cellid = bb.compute_first_entity_collision(point)
        print("expected =", expected_cellid, "found =", found_cellid)
        assert expected_cellid == found_cellid

if __name__ == '__main__':
    test_manifold_point_search()
