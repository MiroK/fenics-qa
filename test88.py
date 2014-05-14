from dolfin import *
import numpy as np

mesh = CircleMesh(Point(0.0,0.0), 1.0, 0.1)
space = FunctionSpace(mesh, "CG", 2)
X = Function(space)
X.interpolate(Expression("x[0]*x[0]+x[1]*x[1]"))

plot(X, interactive=True)

x_coords = np.linspace(-0.9, 0.9, 10)
y_coords = np.zeros(len(x_coords))
z_coords = np.zeros(len(x_coords))

for i, (x, y) in enumerate(zip(x_coords, y_coords)):
  z_coords[i] = X(x, y)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x_coords, y_coords, z_coords)
plt.show()
