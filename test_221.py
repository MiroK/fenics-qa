from dolfin import *
import numpy as np

mesh = UnitSquareMesh(20, 20)
V = FunctionSpace(mesh, 'CG', 2)
f = Expression('sin(x[0]*x[0])+cos(x[1]*x[1])', degree=2)
v = interpolate(f, V)
array = v.vector().array()

dofmap = V.dofmap()
nvertices = mesh.ufl_cell().num_vertices()
# For each cell these are indices of the cell dofs which correspond to vertices
# = entities of dim 0
indices = [dofmap.tabulate_entity_dofs(0, i)[0] for i in range(nvertices)]

vertex_2_dof = dict()
[vertex_2_dof.update(dict(vd for vd in zip(cell.entities(0),
                                           dofmap.cell_dofs(cell.index())[indices])))
 for cell in cells(mesh)]

# For check: the dof values at vertices match those extracted from the array
vertex_indices, vertex_dofs = map(list, zip(*vertex_2_dof.iteritems()))
X = mesh.coordinates()
X = X[vertex_indices]
x, y = X.transpose()
# Dof values
values0 = np.sin(x**2) + np.cos(y**2)
# Extracted
values = array[vertex_dofs]
# Crunch time
print np.linalg.norm(values0-values, np.inf)
