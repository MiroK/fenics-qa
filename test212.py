import dolfin
import numpy as np

N = 4
mesh = dolfin.UnitSquareMesh(N, N, 'crossed')
scheme = 'CR'  # 'TH' or 'CR'
if scheme == 'TH':
    V = dolfin.VectorFunctionSpace(mesh, "CG", 2)
elif scheme == 'CR':
    V = dolfin.VectorFunctionSpace(mesh, "CR", 1)


# define the boundary
def topleftbotright(x, on_boundary):
    return (np.fabs(x[1] - 1.0) < dolfin.DOLFIN_EPS
            or np.fabs(x[0] - 1.0) < dolfin.DOLFIN_EPS
            or np.fabs(x[1]) < dolfin.DOLFIN_EPS
            or np.fabs(x[0]) < dolfin.DOLFIN_EPS)

noslip = dolfin.Constant((0.0, 0.0))
bc = dolfin.DirichletBC(V, noslip, topleftbotright)

bcdict = bc.get_boundary_values()
bcinds = bcdict.keys()

innerinds = np.setdiff1d(range(V.dim()), bcinds)
bcinds = np.setdiff1d(range(V.dim()), innerinds)

# vector of -1
vvec = np.zeros(V.dim()) - 1
# set inner nodes to 1
vvec[innerinds] = 1
# assign to function
v = dolfin.Function(V)
v.vector().set_local(vvec)

# Evaluate CR at midpoints
mesh.init(1, 2)
midpoints = np.array([np.array([facet.midpoint().x(), facet.midpoint().y()])
                      for facet in dolfin.facets(mesh) if facet.exterior()])
U = np.zeros(len(midpoints))
V = np.zeros_like(U)
for i, mp in enumerate(midpoints):
    U[i], V[i] = v(mp)

import matplotlib.tri as tri
import matplotlib.pyplot as plt

# Data for mesh
mesh_coordinates = mesh.coordinates().reshape((-1, 2))
triangles = np.asarray([cell.entities(0) for cell in dolfin.cells(mesh)])
triangulation = tri.Triangulation(mesh_coordinates[:, 0],
                                  mesh_coordinates[:, 1],
                                  triangles)

plt.figure()
plt.triplot(triangulation, color='b')
plt.quiver(midpoints[:, 0], midpoints[:, 1], U, V)
plt.xlim([-0.2, 1.2])
plt.ylim([-0.2, 1.2])
plt.show()
