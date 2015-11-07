from dolfin import *

mesh = UnitSquareMesh(100, 100)

V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

a = inner(u, v)*dx
A = assemble(a)

bc0 = DirichletBC(V, Constant(0), DomainBoundary())
bcs = [bc0, bc0]

[bc.apply(A) for bc in bcs]
