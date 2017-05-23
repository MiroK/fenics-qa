from dolfin import *
import numpy as np

# Setup
mesh = RectangleMesh(Point(-1, -1), Point(1, 1), 128, 128)
cell_domains = CellFunction('size_t', mesh, 0)

solid = '&&'.join(['((0.25 - TOL < x[0]) && (x[0] < 0.75 + TOL))', 
                   '((0.4 - TOL < x[1])  && (x[1] < 0.75 + TOL))'])
solid = CompiledSubDomain(solid, TOL=DOLFIN_EPS)
# Init so that solid point distance to solid is 0
distance_f = VertexFunction('double', mesh, 1)
solid.mark(distance_f, 0)
# Fluid vertices
fluid_vertex_ids = np.where(distance_f.array() > 0.5)[0]
# Represent solid as its own mesh for distance queries
solid.mark(cell_domains, 1)
solid_mesh = SubMesh(mesh, cell_domains, 1)
tree = solid_mesh.bounding_box_tree()
# Fill
for vertex_id in fluid_vertex_ids:
    vertex = Vertex(mesh, vertex_id)
    _, dist = tree.compute_closest_entity(vertex.point())
    distance_f[vertex] = dist
# Let's also build representation as a CG1 function
V = FunctionSpace(mesh, 'CG', 1)
f = Function(V)
transform = dof_to_vertex_map(V)
data = distance_f.array()[transform]
f.vector().set_local(data)
f.vector().apply('insert')

plot(f)
interactive()
