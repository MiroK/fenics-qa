from dolfin import *
import numpy as np

mesh = UnitSquareMesh(20, 20)

points = [Point(0, 0), Point(0.5, 0.5), Point(1, 1)]
n_points = len(points)

V = FunctionSpace(mesh, 'CG', 1)
Q = VectorFunctionSpace(mesh, 'R', 0, dim=n_points)
M = MixedFunctionSpace([V, Q])

# u, p = TrialFunctions(M)
# v, q = TestFunctions(M)
u, v = TrialFunction(V), TestFunction(V)

sources = [PointSource(V, point) for point in points]

a = inner(grad(u), grad(v))*dx
# a = sum(p[i]*inner(v, sources[i])*dx for i in range(n_points))+\
# a = sum(q[i]*inner(u, sources[i])*dx for i in range(n_points))

L = sources[0]
bc = DirichletBC(M.sub(0), Constant(0), DomainBoundary())

uh = Function(V)
solve(a == L, uh, bc)

plot(uh)
interactive()




