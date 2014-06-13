from dolfin import *
import numpy as np

mesh = RectangleMesh(-1, -1, 1, 1, 100, 100)
V = FunctionSpace(mesh, 'CG', 1)

u = Function(V)
U = u.vector()

dofmap = V.dofmap()
dof_x = dofmap.tabulate_all_coordinates(mesh).reshape((-1, 2))
first_dof, last_dof = dofmap.ownership_range()  # U.local_size()

rank = MPI.process_number()
new_values = np.zeros(last_dof - first_dof)
for i in range(len(new_values)):
    x, y = dof_x[i]
    new_values[i] = abs(x)**(rank+1) + abs(y)**(rank+1)

U.set_local(new_values)
U.apply('insert')

plot(u, title=str(rank))
interactive()
