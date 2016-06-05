from dolfin import *

mesh = UnitSquareMesh(20, 20)

f = Expression('1', mpi_comm=mpi_comm_world())
V = FunctionSpace(mesh, 'CG', 1)

v = interpolate(f, V).vector()

print (v.size(), len(v.get_local()))
