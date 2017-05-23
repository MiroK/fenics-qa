from dolfin import *
from mshr import *

center = Point(5, 5)
rad = 5

outer = Rectangle(Point(0,0), Point(10,10))
inner = Circle(center, rad)
domain = outer - inner
domain.set_subdomain(42, inner)

mesh = generate_mesh(domain, 5)

# Get all bdry facets
marked_facets = FacetFunction('size_t', mesh, 0)
DomainBoundary().mark(marked_facets, 1)  # Now ALL bdry (inner & outer) facets are 1

def outer_facet(x): near(x[0]*(10-x[0]), 0) or near(x[1]*(10-x[1]), 0)

# Using how you defined the inner domain mark the inner bdry facets
inner_facets = (facet for facet in SubsetIterator(marked_facets, 1) if
                not outer_facet(facet.midpoint()))
for facet in inner_facets:
    marked_facets[facet] = 2

plot(marked_facets)
interactive()
