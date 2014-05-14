from dolfin import *

mesh = BoxMesh(-1, -1, -1, 1, 1, 1, 20, 20, 20)
h = mesh.hmin()

origin = Point(0., 0., 0.)
inside_cells = CellFunction('size_t', mesh, 0)
for cell in cells(mesh):
  if cell.midpoint().distance(origin) < 0.5 + h:
    insice_cells[cell] = 1
