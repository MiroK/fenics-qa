from dolfin import *
mesh = UnitCubeMesh(2, 2, 2)
d = mesh.topology().dim()

# There's always cell-vertex connectivity
print Cell(mesh, 0).entities(0)

# Not the other way around
print Vertex(mesh, 0).entities(d)

# You must initialize it
mesh.init(0, d)
print Vertex(mesh, 0).entities(d)
