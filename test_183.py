from dolfin import *

mesh = UnitSquareMesh(40, 40, 'crossed')
bmesh = BoundaryMesh(mesh, 'exterior')
# Suppose now bottom abd left boundary are not of interest
cell_f = CellFunction('size_t', bmesh, 0)
CompiledSubDomain('near(x[1], 0)').mark(cell_f, 1)
bmesh_sub = SubMesh(bmesh, cell_f, 0)
tree = bmesh_sub.bounding_box_tree()

V = FunctionSpace(mesh, 'CG', 1)
v_2_d = vertex_to_dof_map(V)
bdry_distance = Function(V)
values = bdry_distance.vector().array()
for index, vertex in enumerate(vertices(mesh)):
    _, d = tree.compute_closest_entity(vertex.point())
    values[v_2_d[index]] = d
bdry_distance.vector().set_local(values)
bdry_distance.vector().apply('insert')

plot(bdry_distance)
interactive()
