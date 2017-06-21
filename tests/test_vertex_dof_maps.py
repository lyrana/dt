from dolfin import *

mesh = UnitSquareMesh(5,5)

# A 2D vector function space represented by CG1 elements
Q = VectorFunctionSpace(mesh, "CG", 1, dim=2)

# DOF indices index the functions in the CG function space
# This map converts a vertex index to a DOF index
v2d=vertex_to_dof_map(Q)
d2v=dof_to_vertex_map(Q)

print "before reshape, v2d=", v2d
v2d = v2d.reshape((-1, mesh.geometry().dim()))
print "after reshape, v2d=", v2d

print "before processing, d2v=", d2v
d2v = d2v[xrange(0, len(d2v), 2)]/2 # the even values, divided by 2??
print "before processing, d2v=", d2v

# Test
# A function on the mesh vertices
# This function takes the index of a vertex and gives a value on that vertex.
v = VertexFunction("size_t", mesh)

# A function represented by CG1 finite-elements basis functions on the mesh
# This function uses the dofs.
q = Function(Q)

# set the value at vertices 10 and 25
v.array()[10] = 1.0
v.array()[25] = 1.0
# 
# Set the dofs at these vertics.
# Sets both x and y components to 1?
q.vector()[v2d[10]] = 1.0
q.vector()[v2d[25]] = 1.0


plot(v)
plot(q)          # press m to visualize the mesh
interactive()
