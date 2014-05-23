from dolfin import *

domain = Rectangle(-1, -1, 1, 1)
mesh = Mesh(domain, 40)
cell_f = CellFunction('size_t', mesh, 0)

circle = lambda x, on_bondary: x[0]*x[0] + x[1]*x[1] < 0.25
AutoSubDomain(circle).mark(cell_f, 1)

circle_mesh = SubMesh(mesh, cell_f, 1)

circle_boundary_mesh = BoundaryMesh(circle_mesh, 'exterior')

