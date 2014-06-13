from dolfin import *
import numpy as np

#domain = Rectangle(0, 0, 1, 1) - Circle(0.5, 0.5, 0.25)
#mesh = Mesh(domain, 20)
mesh = UnitSquareMesh(10, 10)

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

bdry_dofs.sort()

# Consider parallel
bdry_dofs -= first_dof

# Some function to work with
u = Expression('2 + fabs(x[0]-0.5) + fabs(x[1]-0.5)')
u = interpolate(u, V)
U = u.vector()

# Set function's boundary dofs
values = U.get_local()
values[bdry_dofs] = 0
U.set_local(values)

# Check. Project higher order to 0 to make plots same
if order:
    u = project(u, FunctionSpace(mesh, 'DG', 0))
File('u.pvd') << u
plot(u)
interactive()
