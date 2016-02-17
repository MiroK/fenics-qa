from dolfin import *
import numpy as np

N = 4
mesh = UnitSquareMesh(N, N)
x = mesh.coordinates().reshape((-1, 2))

# Vertex with cordinates x[i] has a value in the image at ii[i], jj[i]
h = 1./N
ii, jj = x[:, 0]/h, x[:, 1]/h
ii = np.array(ii, dtype=int)
jj = np.array(jj, dtype=int)

print ii

n = int(sqrt(mesh.num_vertices()))
X, Y = np.meshgrid(np.linspace(0, 1, n), np.linspace(0, 1, n))
image = X**2 + Y**2

image_values = image[ii, jj]

V = FunctionSpace(mesh, 'CG', 1)
image_f = Function(V)
v2d = vertex_to_dof_map(V)
image_f.vector()[:] = image_values[v2d]


import matplotlib.pyplot as plt
plt.figure()
plt.pcolor(image)
plt.show()


V = FunctionSpace(mesh, 'DG', 0)
image_f = interpolate(image_f, V)

plot(image_f)
interactive()
