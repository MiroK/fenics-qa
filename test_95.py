from dolfin import *
mesh = UnitSquareMesh(2,2)
V = FunctionSpace(mesh, "CG", 1)
f = Function(V)
e = Function(V)

f.interpolate(e)
f.interpolate(FunctionAXPY(e, 2.))
