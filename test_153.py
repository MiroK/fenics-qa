from dolfin import *

mesh = UnitSquareMesh( 16, 16 )
deg = 3

V = FunctionSpace(mesh, "CG", deg)
f = Function(V)
f.vector()[:] = 10.

W = VectorFunctionSpace(mesh, 'CG', deg)
u = TrialFunction(W)
v = TestFunction(W)
a = inner(u, v)*dx
L = inner(-grad(f), v)*dx
A, b = assemble_system(a, L)

w = Function(W)
solver = LUSolver()
solver.solve(A, w.vector(), b)

print w.vector()[:].min(), w.vector()[:].max()
