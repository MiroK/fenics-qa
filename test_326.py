from dolfin import *

mesh = UnitSquareMesh(32, 32)

V = FunctionSpace(mesh, 'CG', 1)
W = VectorFunctionSpace(mesh, 'CG', 1)

v = TestFunction(V)
f = Expression(('x[0]-2*x[1]', '2*x[1]'), degree=1), V)

L = assemble(
