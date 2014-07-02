from dolfin import *
mesh = UnitSquareMesh(58, 58)
V = FunctionSpace(mesh, 'CG', 1)
u = project(Expression('x[0]*x[1]'), V)

f = HDF5File(mesh.mpi_comm(), 'u.h5', 'w')
#f = HDF5File('u.h5', 'w')

f.write(mesh, 'mesh')
f.write(u, 'u')
