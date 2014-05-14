from dolfin import *
import numpy as np
mesh = CircleMesh(Point(0.0,0.0),1.0,0.1)
space = VectorFunctionSpace(mesh,"CG",2,dim=3)
X = Function(space)
X.interpolate(Expression(("x[0]","x[1]","x[0]*x[0]+x[1]*x[1]")))
x_coords = []
y_coords = []
z_coords = []
for x in np.linspace(-0.9,0.9,10):
    p=Point(x,0.0)

    a, b, c = X(p)
    print a, b, c

    x_coords.append(a)
    y_coords.append(b)
    z_coords.append(c)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x_coords,y_coords,z_coords)
plt.ylim([-1, 1])
plt.show()
