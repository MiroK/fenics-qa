from dolfin import *
import numpy as np

mesh = UnitSquareMesh(20, 20)
V = FunctionSpace(mesh, 'CG', 1)
u = interpolate(Expression('x[0]*x[1]'), V)

dofmap = V.dofmap()
dof_x = dofmap.tabulate_all_coordinates(mesh).reshape((V.dim(), -1))
x = dof_x[:, 0]
# Extract dof indices where some condition on the first coordinate of dof
# is met
indices = np.where(np.logical_and(x > 0.2, x < 0.4))[0]
# Get coordinates of dof
xs = dof_x[indices]
# Get value of dof
vals = u.vector()[indices]

for x, v in zip(xs, vals):
    print x, v

