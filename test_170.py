from dolfin import *

# solve problem \nabla p = f

mesh = UnitSquareMesh(32, 32)
V = VectorFunctionSpace(mesh, "Lagrange", 1)
P = FunctionSpace(mesh, "Lagrange", 1)

p = TrialFunction(P)
v = TestFunction(V)

a = inner(grad(p), v)*dx

A = assemble(a)
print 'A is %d x %d' % (A.size(0), A.size(1))
