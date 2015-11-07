from dolfin import * 
import numpy as np

mesh = UnitSquareMesh(4, 4)

interface = AutoSubDomain(lambda x, on_boundary: near(x[0], 0.5))
facet_f = FacetFunction('size_t', mesh, 0)
interface.mark(facet_f, 1)
iface_facets = SubsetIterator(facet_f, 1)

mesh.init(1, 2)
iface_cells = set(np.hstack(iface_facet.entities(2)
                            for iface_facet in iface_facets).tolist())

    # First compute connectivity edge <--> vertex for edge, vertex on edge_f
    mesh.init(1, 0)
    topology = mesh.topology()
    # Connect between all edge->c of mesh
    all_edge_vertex_c = topology(1, 0)   
V = FunctionSpace(mesh, 'CG', 1)
dofmap = V.dofmap()

iface_dofs = []
for cell in iface_cells:
    cell_dofs = dofmap.cell_dofs(cell)
    for i, facet in enumerate(facets(Cell(mesh, cell))):
        if facet_f[facet] == 1:
            iface_dofs.extend(cell_dofs[dofmap.tabulate_facet_dofs(i)])

print iface_dofs

