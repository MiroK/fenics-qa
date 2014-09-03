from dolfin import *
import numpy as np

mesh = UnitSquareMesh(250, 250)

VV = VectorFunctionSpace(mesh, 'CG', 1)
V = FunctionSpace(mesh, 'CG', 1)

f = interpolate(Expression(('sin(pi*x[0])', 'cos(pi*x[1])')), VV)
f_norm =\
    interpolate(Expression('sqrt(pow(sin(pi*x[0]), 2) + pow(cos(pi*x[1]), 2))'), V)

timer = Timer('s')
timer.start()

f_x = Function(V)
f_y = Function(V)

f.vector()[:] *= f.vector()

# Split so that f_x = f_x**2, f_y = f_y**2
assigner_VV_to_V = FunctionAssigner([V, V], VV)
assigner_VV_to_V.assign([f_x, f_y], f)

# Split so that f_x = f_x**2, f_y = f_y**2

# f_x will hold |f|**2 = f_x**2 + f_y**2
f_x.vector().axpy(1, f_y.vector())

# f_x will hold |f|
f_x.vector().set_local(np.sqrt(f_x.vector().get_local()))
f_x.vector().apply('')

# Error
f_x.vector().axpy(-1, f_norm.vector())
e = f_x.vector().norm('linf')

time = timer.stop()
if MPI.rank(mpi_comm_world()) == 0:
    print 'Error', e
    print time
