from dolfin import *

mesh = UnitSquareMesh(10, 10)
V = VectorFunctionSpace(mesh, 'CG', 1, 5)
u = TrialFunction(V)
v = TestFunction(V)

C = Identity(5)

a = inner(dot(C, u.dx(0)), v)*dx

