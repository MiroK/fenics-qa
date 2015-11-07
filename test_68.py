from dolfin import *
import numpy as np

mesh = UnitSquareMesh(10, 10)
cell_f = CellFunction('size_t', mesh, 0)
# Define sudomain 1. The other is 0
subd1 = AutoSubDomain(lambda x, on_boundary: x[0] < 0.5+DOLFIN_EPS)
subd1.mark(cell_f, 1)
# Mark the cells
submesh0 = SubMesh(mesh, cell_f, 0)
submesh1 = SubMesh(mesh, cell_f, 1)

V, V1, V0 = map(lambda domain: FunctionSpace(domain, 'CG', 2),
                (mesh, submesh0, submesh1))


def get_restriction_matrix(V, W):
    cell_map = W.mesh().data().array('parent_cell_indices', 2)
	
    dofmapV = V.dofmap()
    dofs_xV = dofmapV.tabulate_all_coordinates(V.mesh()).reshape((-1, 2))
    dofmapW = W.dofmap()
    dofs_xW = dofmapW.tabulate_all_coordinates(W.mesh()).reshape((-1, 2))

    dof_map = {}
    for cell in cells(W.mesh()):
	W_dofs = dofmapW.cell_dofs(cell.index())
	
	V_cell = cell_map[cell.index()]
	V_dofs = dofmapV.cell_dofs(V_cell)

	if np.allclose(dofs_xW[W_dofs], dofs_xV[V_dofs]):
       	    dof_map.update({(k, v) for k, v in zip(W_dofs, V_dofs)}) 

get_restriction_matrix(V, V1)

