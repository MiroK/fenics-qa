from dolfin import *
from mshr import *

L = 2
domain = Cylinder(Point(0, 0, 0), Point(0, 0, L), 1, 1, 20)
mesh = generate_mesh(domain, 30)
# Including top and bottom
mesh = BoundaryMesh(mesh, 'exterior')

cell_f = CellFunction('size_t', mesh, 0)
CompiledSubDomain('!(near(x[2], 0) || near(x[2], L))', L=L).mark(cell_f, 1)

mesh = SubMesh(mesh, cell_f, 1)

plot(mesh)
interactive()
