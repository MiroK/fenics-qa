from dolfin import *

mesh = UnitSquareMesh(3, 3)
V = FunctionSpace(mesh, 'CG', 1)

# Find cells that collide with the point
p = Point(0.42, 0.1234)
tree = mesh.bounding_box_tree()
cells_p = tree.compute_collisions(p)

# Among dofs of those cells find one whose coordinate is closest to p
min_dist, dof = 2, None
dofmap = V.dofmap()
for cell in cells_p:
   dofs = dofmap.cell_dofs(cell)
   dofs_x = dofmap.tabulate_coordinates(Cell(mesh, cell))
   distances = [p.distance(q) for q in map(Point, dofs_x)]
   dist = min(distances)
   if dist < min_dist:
       min_dist = dist
       dof = dofs[distances.index(min_dist)]

# Check that the dof is found correctly
f = interpolate(Expression('sqrt(pow(x[0]-x0, 2)+pow(x[1]-x1, 2))',
                            x0=p[0], x1=p[1]),
                V).vector().array()
print 'Min @ dof', dof, 'with error', abs(min_dist-f[dof])
