from dolfin import *

mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, 'CG', 1)
v = Function(V)

n = FacetNormal(mesh)
h = FacetArea(mesh)
R = avg(h)*(jump(grad(v), n))**2*dS
assemble(R)
