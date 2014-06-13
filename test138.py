from dolfin import *

mesh = UnitSquareMesh(2, 2)

# Init facet-cell connectivity
tdim = mesh.topology().dim()
mesh.init(tdim - 1, tdim)

# For every cell, build a list of cells that are connected to its facets
# but are not the iterated cell
cell_neighbors = {cell.index(): sum((filter(lambda ci: ci != cell.index(),
                                            facet.entities(tdim))
                                    for facet in facets(cell)), [])
                for cell in cells(mesh)}

print cell_neighbors
