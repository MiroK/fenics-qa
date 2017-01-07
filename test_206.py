from dolfin import *

r = 0.5
mesh_0 = RectangleMesh(Point(-1, -1), Point(2, 2), 16, 16)
mesh_1 = RectangleMesh(Point(-2, -2), Point(1, 1), 8, 8)

# Build multimesh
multimesh = MultiMesh()
multimesh.add(mesh_0)
multimesh.add(mesh_1)
multimesh.build()

V = MultiMeshFunctionSpace(multimesh, 'DG', 0)
f = MultiMeshFunction(V)
f.vector()[:] = 1.
# compute integral
L = f*dO
vol = assemble_multimesh(L)

print vol
