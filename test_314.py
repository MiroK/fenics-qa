from dolfin import *
from mshr import *

size = 20

# Create mesh 
domain = Rectangle(Point(0., 0.), Point(pi, pi))
mesh = generate_mesh(domain, 20)


class Horizontal(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and ((abs(x[1] - 0.) < DOLFIN_EPS) or (abs(x[1] - pi) < DOLFIN_EPS))

class Vertical(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and ((abs(x[0] - 0.) < DOLFIN_EPS) or (abs(x[0] - pi) < DOLFIN_EPS))


boundaries = FacetFunction("size_t", mesh, 0) 

Horizontal().mark(boundaries, 1) 
Vertical().mark(boundaries, 2) 

V = FunctionSpace(mesh, 'CG', 1)
for i in (1, 2):
    count = sum(1 for _ in SubsetIterator(boundaries, i))
    print 'Facet marked as %d is %d' % (i, count)
    assert len(DirichletBC(V, Constant(0.), boundaries, i).get_boundary_values()) > 0


