import matplotlib.pyplot as plt
import matplotlib.tri as tri
from dolfin import *
import numpy as np

domain = Rectangle(-1, -1, 1, 1) - Circle(0, 0, 0.5)
mesh = Mesh(domain, 20)
n = mesh.num_vertices()
d = mesh.geometry().dim()

# Create the triangulation
mesh_coordinates = mesh.coordinates().reshape((n, d))
triangles = np.asarray([cell.entities(0) for cell in cells(mesh)])
triangulation = tri.Triangulation(mesh_coordinates[:, 0],
                                  mesh_coordinates[:, 1],
                                  triangles)

# Plot the mesh
plt.figure()
plt.triplot(triangulation)
plt.savefig('mesh.svg')

# Create some function
V = FunctionSpace(mesh, 'CG', 2)
f_exp = Expression('sin(2*pi*(x[0]*x[0]+x[1]*x[1]))')
f = interpolate(f_exp, V)

cmap = plt.cm.rainbow
# Get the z values as face colors for each triangle(midpoint)
plt.figure()
zfaces = np.asarray([f(cell.midpoint()) for cell in cells(mesh)])
plt.tripcolor(triangulation, facecolors=zfaces, edgecolors='k', cmap=cmap)
plt.savefig('f0.svg')

# Get the z values for each vertex
plt.figure()
z = np.asarray([f(point) for point in mesh_coordinates])
plt.tripcolor(triangulation, z, cmap=cmap)  # alt plt.tricontourf(...)
plt.colorbar()
plt.savefig('f1.svg')

# Plot of vector field
W = VectorFunctionSpace(mesh, 'DG', 1)
gradf = project(grad(f), W)
z = np.asarray([gradf(point) for point in mesh_coordinates])

z_mag = np.sqrt(np.sum((z*z), axis=1))
plt.figure()
plt.quiver(mesh_coordinates[:, 0], mesh_coordinates[:, 1], z[:, 0], z[:, 1],
           z_mag, cmap=cmap)
plt.colorbar()
plt.savefig('gradf.svg')

# Comment to avoid pop-up
plt.show()
