from dolfin import *
import numpy as np

mesh = RectangleMesh(Point(-1, -1), Point(1, 1), 2, 10)
iface = FacetFunction('size_t', mesh, 0)
CompiledSubDomain('near(x[0], 0)').mark(iface, 1)

mesh.init(1, 0)      
ivertex = list(set(sum((facet.entities(0).tolist()
                        for facet in SubsetIterator(iface, 1)), [])))

dx = np.linspace(0, 2, 10)
dx = 0.5*np.sin(2*pi*dx)

coordinates = mesh.coordinates()
for dxi in dx:
    coordinates[ivertex, 0] = dxi
    plot(mesh, interactive=True)



jvertex = []
# For each interface facet get the vertices connected to it
for facet in SubsetIterator(iface, 1):
    jvertex.extend(facet.entities(0))
# From the construction there are duplicate vertices in the list. Filter them by
# creating a set from the list
jvertex = set(jvertex)
# Bacause we want to use jvertex for indexing arrays it needs to be converted
# back to list
jvertex = list(jvertex)

print jvertex == ivertex
