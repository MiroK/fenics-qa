from scipy.interpolate import LinearNDInterpolator
from dolfin import *
from mshr import *
import numpy as np

# Grid
x, y = np.mgrid[0:2:0.1, 0:5:0.1]
# Pretend data
z = np.sin(x+y)
# Spline interpolant
interpolant = LinearNDInterpolator(np.c_[x.flatten(), y.flatten()], z.flatten())
# Interpolate to FE-space. Unstrucuted mesh
domain = Rectangle(Point(0., 0.), Point(2., 5.))
mesh = generate_mesh(domain, 4)

V = FunctionSpace(mesh, 'CG', 1)
XY = V.dofmap().tabulate_all_coordinates(mesh).reshape((-1, 2))

g = Function(V)
# Expansion coefs of nodal interpolant by eval spline
g_values = interpolant(XY)
print XY
print '>>', interpolant([[2, 0.125],
                         [0, 0.125],
                         [2., 0.125]])

g.vector()[:] = g_values

# Some check
f = interpolate(Expression('sin(x[0]+x[1])'), V)
f.vector().axpy(-1, g.vector())

for i, val in enumerate(g.vector().array()):
    print XY[i], val
