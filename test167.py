from dolfin import *

P_0 = Point(0.9, 0.1)
P_1 = Point(0.1, 0.9)
P_out = Point(10, 10)

mesh = UnitSquareMesh(2, 2)
mesh_f = CellFunction('size_t', mesh, 0)
for cell in cells(mesh):
    M = cell.midpoint()
    x, y = M.x(), M.y()
    if y > x:
        mesh_f[cell] = 1

mesh0 = SubMesh(mesh, mesh_f, 0)
mesh1 = SubMesh(mesh, mesh_f, 1)

tree0 = mesh0.bounding_box_tree()
tree0.build(mesh0, 2)

tree1 = mesh1.bounding_box_tree()
tree1.build(mesh1, 2)

print tree0.compute_first_collision(P_0), tree1.compute_first_collision(P_0)
print tree0.compute_first_collision(P_1), tree1.compute_first_collision(P_1)
print tree0.compute_first_collision(P_out), tree1.compute_first_collision(P_out)
