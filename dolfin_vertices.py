# from dolfin import *
# 
# mesh = Mesh('dolfin_fine.xml.gz')
# # We will be getting vertices from facets that are on the boundary of domain
# facet_f = FacetFunction('size_t', mesh, 0)
# 
# # Mark all facets on the boundary as 1
# DomainBoundary().mark(facet_f, 1)
# 
# on_outer_boundary = lambda x: near(x[0]*(1-x[0]), 0) or near(x[1]*(1-x[1]), 0)
# # Mark facets on the boundary of the square as 2
# AutoSubDomain(lambda x, on_boundary:\
#               on_boundary and on_outer_boundary(x)).mark(facet_f, 2)
# 
# # Get the vertices connected to facets marked with 2
# mesh.init(1, 0)
# dolfin_vertices = set(sum((facet.entities(0).tolist()
#                            for facet in SubsetIterator(facet_f, 2)), []))
# 
# # Check
# vertex_x = mesh.coordinates().reshape((-1, 2))
# assert all(on_outer_boundary(vertex_x[vertex]) for vertex in dolfin_vertices)

from dolfin import Point, plot, interactive
from mshr import *

domain = Circle(Point(0, 0), 1) - Circle(Point(0.5, 0.5), 0.2)\
                                - Circle(Point(-0.5, 0.5), 0.2)\
                                - Circle(Point(0, -0.5), 0.35)\
                                - Rectangle(Point(-0.125, -0.125), Point(0.125, 0.5))

mesh = generate_mesh(domain, 50)

plot(mesh, title='What!')
interactive()



