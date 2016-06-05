from dolfin import *
import sys

n_cells = int(sys.argv[1])

mesh = UnitSquareMesh(n_cells, n_cells)                                                   
V = FunctionSpace(mesh, 'CG', 2)  
u = Function(V)                                                                 

vec = u.vector()
values = vec.get_local()                                            

dofmap = V.dofmap()                                                             
my_first, my_last = dofmap.ownership_range()                # global

# 'Handle' API change of tabulate coordinates
if dolfin_version().split('.')[1] == '7':
    x = V.tabulate_dof_coordinates().reshape((-1, 2))
else:
    x = V.dofmap().tabulate_all_coordinates(mesh)

unowned = dofmap.local_to_global_unowned()
dofs = filter(lambda dof: dofmap.local_to_global_index(dof) not in unowned, 
              xrange(my_last-my_first))
x = x[dofs]

import numpy as np
values[:] = x[:, 0]**2 + x[:, 1]**3

vec.set_local(values)
vec.apply('insert')

# Check
u0 = interpolate(Expression('x[0]*x[0]+x[1]*x[1]*x[1]', degree=2), V)
u0.vector().axpy(-1, vec)

error = u0.vector().norm('linf')
if MPI.rank(mpi_comm_world()) == 0: print V.dim(), error
