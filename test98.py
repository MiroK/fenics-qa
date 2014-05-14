from dolfin import *

mesh = UnitSquareMesh(50, 50)
V = FunctionSpace(mesh, 'CG', 1)

u = TrialFunction(V)
v = TestFunction(V)

f = Expression('sin(2*pi*x[0])*sin(pi*x[1])')
a = inner(grad(u), grad(v))*dx
L = inner(f, v)*dx
bc = DirichletBC(V, Constant(0), DomainBoundary())

phi = Function(V)
solve(a == L, phi, bc)

grad_norm = norm(phi, 'h10')
psi0 = project(phi/grad_norm, V, bcs=bc)

plot(psi0, interactive=True)

