from dolfin import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri

mesh = UnitSquareMesh(32, 32)
V = FunctionSpace(mesh, 'CG', 1)
v = interpolate(Expression('std::abs(x[0]+x[1])', degree=2), V)

triang = tri.Triangulation(*mesh.coordinates().reshape((-1, 2)).T,
                           triangles=mesh.cells())
Z = v.compute_vertex_values(mesh)

plt.figure()
plt.tricontourf(triang, Z)
plt.colorbar()
plt.show()
