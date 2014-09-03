from dolin import *

mesh = UnitIntervalMesh(10)
V = FunctionSpace(mesh, 'CG', 1)

u = TrialFunction(V)
v = TestFunction(V)
a = inner(grad(u), grad(v))*dx


