from dolfin import *
from random import sample

mesh = UnitIntervalMesh(20)
n_vertices = mesh.num_vertices()
# Get all vertices in random order
data = sample(range(n_vertices), n_vertices)
# The values are such that we make function f=x
values = mesh.coordinates().reshape((n_vertices, ))[data]

V = FunctionSpace(mesh, 'CG', 1)
f = Function(V).vector()

# First approach.
v2d = vertex_to_dof_map(V)
for vertex, value in zip(data, values):
    dof = v2d[vertex]
    f[dof] = value

# Compare to f0=x
f0 = interpolate(Expression('x[0]'), V).vector()
f0.axpy(-1, f)
print f0.norm('linf')

# Second approach. We will use it to zero the vector by prescribing 
# negative values
gdim = mesh.geometry().dim()
for vertex, value in zip(data, values):
    v = Vertex(mesh, vertex)
    x = Point(*[v.x(i) for i in range(gdim)])
    p = PointSource(V, x, -value)
    p.apply(f)

# Compare to 0
print f.norm('linf')
