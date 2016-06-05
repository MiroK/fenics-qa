from dolfin import *
import numpy as np

mesh = UnitSquareMesh(4, 4)
f = EdgeFunction('size_t', mesh, 0)
DomainBoundary().mark(f, 1)

# Get ids for 1 marked edges
before = [edge.index() for edge in SubsetIterator(f, 1)]
# and physical coordinates of their midpoints
before_x = [Edge(mesh, edge).midpoint() for edge in before]
# Deform
ALE.move(mesh, Constant((1, 1)))
# The mesh has moved but ids have not changed
after = set(np.where(f.array() == 1)[0])
assert len(before) == len(after)
assert all(edge in after for edge in before)
# Only the physical coordinates of the midpoints did
after_x = [Edge(mesh, edge).midpoint() for edge in before]
print min(a.distance(b) for a, b in zip(after_x, before_x))
