from mshr import Sphere, generate_mesh
from dolfin import *
import numpy as np

# Radius of outer and inner sphere
oradius, iradius = 5., 1.

# Geometry
outer_sphere = Sphere(Point(0., 0., 0.), oradius)
inner_sphere = Sphere(Point(0., 0., 0.), iradius)
g3d = outer_sphere - inner_sphere
mesh = generate_mesh(g3d, 8)

facet_f = FacetFunction('size_t', mesh, 0)
# Mark all external facets
DomainBoundary().mark(facet_f, 1)
# Now mark the inner boundary by checking that radius is smaller than radius
# of mid sphere
center = Point(0., 0., 0.)
mid = 0.5*(oradius - iradius)
for facet in SubsetIterator(facet_f, 1):
    if facet.midpoint().distance(center) < mid:
        facet_f[facet] = 2

# See that both values are now present
print len(np.where(facet_f.array() == 1)[0])
print len(np.where(facet_f.array() == 2)[0])
