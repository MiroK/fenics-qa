from dolfin import *
from mshr import *

def foo(N):
    P = Point(0, 0)
    circle = Circle(P, 1)
    mesh = generate_mesh(circle, N)
    mesh = BoundaryMesh(mesh, 'exterior')
    return mesh.num_cells()

import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(5, 135, 5)
y = map(foo, x)

A = np.vstack([x, np.ones_like(x)]).T
m, c = np.linalg.lstsq(A, y)[0]

#plt.figure()
#plt.plot(x, y)
#plt.plot(x, m*x+c, 'r--')
#plt.show()

N = 100
P = Point(0, 0)
circle = Circle(P, 1)
mesh = generate_mesh(circle, N)
mesh = BoundaryMesh(mesh, 'exterior')
h = mesh.hmin()

h_predict = 2*pi/(m*N+c)

print h, h_predict
