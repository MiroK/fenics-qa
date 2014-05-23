from dolfin import *

mesh = UnitSquareMesh(20, 20)
V = FunctionSpace(mesh, 'CG', 4)
dm = V.dofmap()

f = Expression('sin(pi*x[0])*cos(pi*x[1])')
f = interpolate(f, V)
F = f.vector()

# Use edge function to store values
F_min = EdgeFunction('double', mesh, F.max())
F_max = EdgeFunction('double', mesh, F.min())

for cell in cells(mesh):
    cell_index = cell.index()
    cell_dofs = dm.cell_dofs(cell_index)

    for i, edge in enumerate(edges(cell)):
        # These are dofs internal to edge
        edge_dofs = list(cell_dofs[dm.tabulate_facet_dofs(i)])

        # Add vertex dofs
        vertex_dofs = [cell_dofs[dm.tabulate_facet_dofs(cell.index(v))][0]
                       for v in vertices(edge)]

        # Combine
        dofs = edge_dofs + vertex_dofs

        # Find min/max
        F_min[edge] = min(F[dofs].min(), F_min[edge])
        F_max[edge] = max(F[dofs].max(), F_max[edge])

plot(f)
plot(F_min, title='Min')
plot(F_max, title='Max')
interactive()
