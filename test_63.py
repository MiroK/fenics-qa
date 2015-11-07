from dolfin import *
from mshr import *

domain = Circle(Point(0, 0), 1.0)
domain.set_subdomain(1, Circle(Point(0, 0.8), 0.05))
mesh = generate_mesh(domain, 40)
cell_f = MeshFunction("size_t", mesh, 2, mesh.domains())

plot(cell_f, interactive=True)
