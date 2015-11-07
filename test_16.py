from dolfin import *
from mshr import *
import numpy as np

res = 10
mesh_type = 'mshr'

get_mesh = lambda domain, res:\
           generate_mesh(domain, res) if mesh_type == 'mshr' else\
          UnitSquareMesh(20, 20)

domain = Rectangle(Point(0, 0), Point(1, 1)) - Circle(Point(0.5, 0.5), 0.125)

x_ = get_mesh(domain,res).coordinates().reshape((-1, 2))

n_mesh_changes = 0
for i in range(100):
    x = get_mesh(domain,res).coordinates().reshape((-1, 2))
    if np.linalg.norm(x - x_) > 1E-14:
        n_mesh_changes += 1
    np.copyto(x_, x)

print 'Mesh changed', n_mesh_changes

class PeriodicBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS
    
    def map(self, x, y):
        y[0] = x[0] - 1.0
        y[1] = x[1]

pbc = PeriodicBoundaryComputation()
mesh = get_mesh(domain,res)
mesh.init(1)

paired_edges = pbc.compute_periodic_pairs(mesh, PeriodicBoundary(), 1)
edge_f = FacetFunction('size_t', mesh, 0)

for i, (edge0, edge1) in enumerate(paired_edges.iteritems(), 1):
    edge_f[edge0] = i
    edge_f[edge1[1]] = i

plot(edge_f)
interactive()


