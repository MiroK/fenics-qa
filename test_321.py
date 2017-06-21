from dolfin import *

mesh = UnitSquareMesh(32, 32)
V = VectorFunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)
f = Constant((1, 1))

a = inner(u, v)*ds
L = inner(f, v)*ds
A = assemble(a, keep_diagonal=True)
b = assemble(L)
A.ident_zeros()

uh = Function(V)
solve(A, uh.vector(), b)

plot(uh, interactive=True)
