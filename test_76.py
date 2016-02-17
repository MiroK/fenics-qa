from dolfin import *
import numpy as np

mesh = UnitSquareMesh(3, 3)
V = FunctionSpace(mesh, 'CG', 3)
f = interpolate(Expression('3*x[0]*x[0]*x[1]-2*x[1]*x[1]*x[0]'), V)
el = V.element()

# Where to evaluate
x = np.array([0.33, 0.55])

# Find the cell with point
x_point = Point(*x) 
cell_id = mesh.bounding_box_tree().compute_first_entity_collision(x_point)
cell = Cell(mesh, cell_id)
coordinate_dofs = cell.get_vertex_coordinates()

# Array for values = single scalar here
values = np.zeros(1, dtype=float)
# Compute all 2nd order derivatives
for i in range(el.space_dimension()):
    el.evaluate_basis(i, values, x, coordinate_dofs, cell.orientation())
    print i, values
