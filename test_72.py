from dolfin import *

mesh = UnitSquareMesh(4, 4)
V = FunctionSpace(mesh, 'CG', 1)
W = MixedFunctionSpace([V, V])

v1 = Function(V)
v2 = Function(V)
w = Function(W)
w.vector()[:] = 1

fa = FunctionAssigner([V, V], W)
fa.assign([v1, v2], w)

print v2.vector().array()
print v2.vector().copy()
