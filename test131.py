from dolfin import *

mesh = UnitSquareMesh(100, 100)
V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

f = Expression('sin(pi*x[0])*cos(pi*x[1])')
a = inner(grad(u), grad(v))*dx
L = inner(f, v)*dx
bc = DirichletBC(V, Constant(0), DomainBoundary())

A, b = assemble_system(a, L, bc)
uh = Function(V)

solve(A, uh.vector(), b, 'gmres', 'hypre_euclid')

plot(uh)
interactive()
