from dolfin import *
from numpy import zeros

mesh = UnitSquareMesh(10, 10)

V = FunctionSpace(mesh, 'CG', 1)
u = Function(V)

plot(u, interactive=True)

for i in range(V.dim()):
    values = zeros(V.dim())
    values[i] = 1
    u.vector()[:] = values
    plot(u, interactive=True)
