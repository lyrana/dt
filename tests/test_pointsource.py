__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

from dolfin import *

class Delta(Expression):
    def __init__(self, eps):
        self.eps = eps
    def eval(self, values, x):
        eps = self.eps
        values[0] = eps/pi/(x[0]**2 + eps**2) 

mesh = IntervalMesh(1000, -1, 1)
V = FunctionSpace(mesh, 'CG', 1)
bc = DirichletBC(V, Constant(0), 'on_boundary')
u = TrialFunction(V)
v = TestFunction(V)

a = inner(grad(u), grad(v))*dx 
delta = Delta(1E-4)
L = inner(delta, v)*dx

u = Function(V)
if False:
    solve(a == L, u, bc)
else:
    A, b = assemble_system(a, L, bc)
    b.zero()
    PointSource(V, Point(0), 100).apply(b)
    solve(A, u.vector(), b)

plot(delta, mesh=mesh)
plot(u, interactive=True)
