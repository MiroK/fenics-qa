from dolfin import *
import numpy as np
import sys

field = sys.argv[1]

mesh = UnitCubeMesh(2, 2, 2)
bmesh = BoundaryMesh(mesh, 'exterior')

if field == 'scalar':
    V = FunctionSpace(mesh, 'CG', 1)
    Vb = FunctionSpace(bmesh, 'CG', 1)
    F = Expression('x[0]+x[1]', degree=1)
else:
    V = VectorFunctionSpace(mesh, 'CG', 1)
    Vb = VectorFunctionSpace(bmesh, 'CG', 1)
    F = Expression(('2*(x[0]+x[1]) + 4*(x[1]-x[2]) - x[2]',
                    'x[0]+x[1]',
                    'x[1]+x[2]'), degree=1)

dofmapb = Vb.dofmap()
firstb, lastb = dofmapb.ownership_range()
dofmap = V.dofmap()
first, last = dofmap.ownership_range()

Xb = Vb.tabulate_dof_coordinates().reshape((-1, 3))
X = V.tabulate_dof_coordinates().reshape((-1, 3))

cb2f = bmesh.entity_map(2)
mesh.init(2, 3)
mesh.init(3, 2)
f2c = mesh.topology()(2, 3)
c2f = mesh.topology()(3, 2)
for i, cellb in enumerate(cells(bmesh)):
    dofsb = dofmapb.cell_dofs(i)

    facet = cb2f[i]                                 # Global
    cell = int(f2c(facet))  
    facet = c2f(cell).tolist().index(facet)         # Local
    dofs = dofmap.cell_dofs(cell)                   # All
    dofs = dofs[dofmap.tabulate_facet_dofs(facet)]  # On the facet

    # Make sure that they match
    info('%s %s' % ([firstb <= dofmapb.local_to_global_index(d) < lastb for d in dofsb],
                    [first <= dofmap.local_to_global_index(d) < last for d in dofs]))
    #xb = Xb[dofsb]
    #x = X[dofs]
    #print [np.linalg.norm(xi - xbi) for xi, xbi in zip(x, xb)]
    print Facet(mesh, facet).get_vertex_coordinates() - cell_b.get_vertex_coordinates()

    
