from dolfin import *
import numpy as np

mesh = UnitSquareMesh(64, 64)
facet_f = FacetFunction('size_t', mesh, 0)
CompiledSubDomain('near(x[0], x[1])').mark(facet_f, 1)

V = FunctionSpace(mesh, 'CG', 1)
v = interpolate(Expression('x[0]+x[1]', degree=1), V)

mesh.init(1, 0)
# Unique vertices of the line (by indices)
vs = list(set(sum((f.entities(0).tolist() for f in SubsetIterator(facet_f, 1)), [])))
# Get degrees of freedom associated to vertices (by indices)
v2d = vertex_to_dof_map(V)
d = v2d[vs]
# Dof values
values = v.vector().array()[d]
# The values should really be the same as evaluating v at vertices (by coordinates)
x = mesh.coordinates().reshape((-1, 2))
# Collect values at points that are on line
values0 = np.array([v(xi) for xi in x[vs]])

print np.linalg.norm(values0 - values)
