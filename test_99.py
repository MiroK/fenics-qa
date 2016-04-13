import dolfin
from mshr import *

dolfin.set_log_level(dolfin.TRACE)

# Define 2D geometry
domain = Rectangle(dolfin.Point(0., 0.), dolfin.Point(1., 1.)) - Circle(dolfin.Point(0.0, 0.0), .35)
# domain.set_subdomain(1, Rectangle(dolfin.Point(.05, .05), dolfin.Point(.95, .95)))
domain.set_subdomain(2, Circle(dolfin.Point(0, 0), .45))
# domain.set_subdomain(3, Circle(dolfin.Point(0,0), .6))

# Generate and plot mesh
mesh2d = generate_mesh(domain, 45)
dolfin.plot(mesh2d, "2D mesh")

# Convert subdomains to mesh function for plotting
# mf = dolfin.MeshFunction("size_t", mesh2d, 2, mesh2d.domains())
# dolfin.plot(mf, "Subdomains")

dolfin.interactive()
