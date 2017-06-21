from dolfin import *
from mshr import *

domain = Sphere(Point(0, 0, 0), 1)
mesh = generate_mesh(domain, 50)

n = FacetNormal(mesh)
print assemble(div(n)*ds(domain=mesh))
