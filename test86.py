from dolfin import *

mesh = RectangleMesh(-1, -1, 1, 1, 100, 100)
V = FunctionSpace(mesh, "CG", 1)
u = TrialFunction(V)
v = TestFunction(V)

a = inner(grad(u), grad(v))*dx
L = Constant(0)*v*dx
bc = DirichletBC(V, Constant(0), DomainBoundary())
A, b = assemble_system(a, L, bc)

delta = PointSource(V, Point(0., 0.,), 100)
delta.apply(b)

u = Function(V)
solve(A, u.vector(), b)

plot(u, interactive=True)
