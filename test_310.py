from dolfin import *
import numpy as np

# Create mesh and define function space
mesh = RectangleMesh(Point(0,0),Point(1,1),10,10)
Mh = FunctionSpace(mesh, "Lagrange", 1)

# Tabulating coordinates
# dof_coord = Mh.dofmap().tabulate_all_coordinates(mesh)   
dof_coord = Mh.tabulate_dof_coordinates().reshape((-1, 2))  # 2016.2.0

# Create a function with maximum near the center
f = Expression("1.0/exp(pow(x[0] - 0.5,2)+pow(x[1]-0.5,2))", degree=3)
v = Function(Mh); v=interpolate(f,Mh)
plot(v, interactive=True)

# Getting the indexes of the values above a certain threshold 
indexes = np.where(v.vector()>0.9)[0]
values  = v.vector().array()[indexes]

# Printing the value and locations 
for j, dof_index in enumerate(indexes):
    x, y = dof_coord[dof_index]
    val = values[j]
    val0 = v(x, y)
    print 'v(%2g,%2g) = %g [%g]' % (x, y, val, val0)


