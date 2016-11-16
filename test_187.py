from dolfin import *

mesh = UnitSquareMesh(10, 10)
bmesh = BoundaryMesh(mesh, 'exterior')
# Map vertices of bmesh to vertices of mesh
vb_2_v = bmesh.entity_map(0)
# Check
assert all(near(Vertex(bmesh, vb_index).point().distance(Vertex(mesh, v_index).point()), 0)
           for vb_index, v_index in enumerate(vb_2_v))
# The reverse map
v_2_vb = {v: vb for vb, v in enumerate(vb_2_v)}

