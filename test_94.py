from dolfin import *
import numpy as np

mesh = UnitSquareMesh(2, 2)
V = FunctionSpace(mesh, 'CG', 1)
first, last = V.dofmap().ownership_range()

comm = mpi_comm_world()
rank = MPI.rank(comm)
bc = DirichletBC(V, Constant(rank), 'on_boundary')
print rank, 'x', bc.get_boundary_values().keys()
dofs = filter(lambda dof: first <= dof < last, bc.get_boundary_values().keys())
print rank, 'x', dofs
dofs = np.array(dofs) - first
print rank, 'x', dofs

u = TrialFunction(V)
v = TestFunction(V)
a = Constant(0)*inner(u, v)*dx
A = assemble(a, keep_diagonal=True)

diag = Vector(comm, A.size(0))
A.init_vector(diag, 0)

diag_values = np.zeros(diag.local_size())
diag_values[dofs] = 2*(rank + 1)
diag.set_local(diag_values)
diag.apply('insert')
A.set_diagonal(diag)

print A.array()
