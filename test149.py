from dolfin import *

mesh = UnitSquareMesh(10, 10)
gdim = mesh.geometry().dim()
V = FunctionSpace(mesh, 'DG', 1)
dofmap = V.dofmap()

dofs = dofmap.dofs()
# Get coordinates as len(dofs) x gdim array
dofs_x = dofmap.tabulate_all_coordinates(mesh).reshape((-1, gdim))

for dof, dof_x in zip(dofs, dofs_x):
	print dof, ':', dof_x

# See how many dofs are associated with vertices 
print dofmap.num_entity_dofs(0)
# All the dofs are associated with cells
print dofmap.num_entity_dofs(2), V.dim()/mesh.num_cells()

