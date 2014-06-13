from dolfin import *

mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, 'CG', 1)

u = Function(V)
f = interpolate(Expression('sin(x[0])'), V)

u.assign(u + f)
