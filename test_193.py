from dolfin import *
import numpy as np

mesh = UnitSquareMesh(4, 4)
V1 = VectorFunctionSpace(mesh, 'CG', 1)
V2 = VectorFunctionSpace(mesh, 'CG', 2)

f = Expression(('x[0]', '-x[1]'), degree=1)
v2 = interpolate(f, V2)
vertex_values = v2.compute_vertex_values()

N = 2

v2d = vertex_to_dof_map(V1)
d2v = dof_to_vertex_map(V1)

v1 = Function(V1)
values = v1.vector().get_local()

first, last = V1.dofmap().ownership_range()
mask = np.array(map(lambda d: first <= V1.dofmap().local_to_global_index(d) < last, v2d))
v2d = v2d[mask]
vertex_value = (vertex_values[d2v].reshape((N, -1)).T.flatten())[mask]

for i in range(N):
    values[v2d[i::N]] = vertex_values[v2d[i::N]]

v1.vector().set_local(values)
v0 = interpolate(f, V1)

plot(v0)
plot(v1)
interactive()
# 
print (v0.vector()-v1.vector()).norm('linf')
