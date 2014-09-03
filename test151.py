from dolfin import *
mesh = UnitSquareMesh(58, 58)
print type(mesh)
V = FunctionSpace(mesh, 'CG', 1)
u = project(Expression('x[0]*x[1]'), V)

def write_solution(ouput_filename):
    g = HDF5File(mesh.mpi_comm(), ouput_filename, 'w')
    g.write(u,'u')

#f.write(mesh, 'mesh')
#f.write(u, 'u')
