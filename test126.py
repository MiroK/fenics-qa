from dolfin import *

mesh = CircleMesh(Point(0, 0), 1, 0.5)

f = Constant(1.)
S = FunctionSpace(mesh, 'CG', 1)

n = FacetNormal(mesh)

V = FunctionSpace(mesh, 'RT', 1)
v = project(f*n, V)
