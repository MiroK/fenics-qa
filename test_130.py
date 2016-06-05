from mshr import *
from dolfin import *

domain = Circle(Point(0, 0), 1)
domain.set_subdomain(

mesh = generate_mesh(domain, 40)
plot(mesh)
interactive()
