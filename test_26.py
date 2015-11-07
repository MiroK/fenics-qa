from fenics import *

mesh = UnitSquareMesh(20, 20)
V = FunctionSpace(mesh, 'CG', 2)

cell_f = CellFunction('size_t', mesh, 0)
AutoSubDomain(lambda x: x[0] > x[1]).mark(cell_f, 1)
domain_cells = SubsetIterator(cell_f, 1)

dofmap = V.dofmap()
# Dofs in domain
dofs = sum((dofmap.cell_dofs(cell.index()).tolist()
           for cell in domain_cells), [])
# Get unique dofs in local numbering that the process sees. Some might not be owned
dofs = set(dofs)
# Where in global vector are dofs owned by the process
my_first, my_last = dofmap.ownership_range()
# Keep only owned ones, ie those whose global index is in ownership range
dofs = filter(lambda dof: my_first <= dofmap.local_to_global_index(dof) < my_last,
              dofs)
# Assign
f = Function(V)
f.vector()[dofs] = 1

plot(f)
interactive()
