from dolfin import *

mesh = UnitSquareMesh(10, 10)

p = 0  # p > 0 okay
V = FunctionSpace(mesh, 'DG', p)

u = interpolate(Expression('x[0]'), V)

g = u.compute_vertex_values(mesh)
#print len(g), mesh.num_vertices()
#print g

plot(u)
interactive()#
