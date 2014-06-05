from dolfin import *
import numpy as np

# 2d test
domain = Rectangle(0, 0, 1, 1) - Circle(0.5, 0.5, 0.25)

mesh = Mesh(domain, 20)

# Get cell facet connectivity
tdim = mesh.topology().dim()
mesh.init(tdim, tdim-1)

# Mark facets on the boundary
bdry_facets = FacetFunction('bool', mesh, False)
BDRY_FACETS = bdry_facets.array()
DomainBoundary().mark(bdry_facets, True)

# Get all dofs which belong to cells with some facet on the boundary
order = 3
V = FunctionSpace(mesh, 'DG', order)
dofmap = V.dofmap()
first_dof, last_dof = dofmap.ownership_range()

bdry_dofs = np.concatenate([dofmap.cell_dofs(cell.index())
                            for cell in cells(mesh)
                            if any(BDRY_FACETS[cell.entities(tdim - 1)])])

bdry_dofs -= first_dof
bdry_dofs.sort()

# Set function's boundary dofs
u = Function(V)
U = u.vector()
values = np.zeros(last_dof - first_dof)
values[bdry_dofs] = 1
U.set_local(values)

# Check,
# Project higher order to 0 to make plots same
if order:
    u = project(u, FunctionSpace(mesh, 'DG', 0))
File('uu.pvd') << u
