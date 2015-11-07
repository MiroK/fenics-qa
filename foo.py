from dolfin import *
import numpy
from matplotlib import pyplot

import logging
logging.getLogger('FFC').setLevel(logging.WARNING)
logging.getLogger('UFL').setLevel(logging.WARNING)

# Make mesh ghosted to permit evaluation of DG terms
# ("shared_facet", "shared_vertex" or "none")
# http://fenicsproject.org/pipermail/fenics-support/2014-July/000778.html
parameters["ghost_mode"] = "shared_vertex"

class PeriodicDomain(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0)

    def map(self, x, y):
        y[0] = x[0] - 1.0

constrained_domain = PeriodicDomain()

mesh = UnitIntervalMesh(3)
V0 = FunctionSpace(mesh, 'DG', 1)
V1 = FunctionSpace(mesh, 'DG', 1, constrained_domain=constrained_domain)

As = []
vel = 0.1
dt = 0.1
for V in (V0, V1):
    print V.dim()


    u = TrialFunction(V)
    v = TestFunction(V)

    uf = Function(V)  # New value
    upf = Function(V) # Previous value

    # u0 = Expression(u0code)
    # project(u0, V, function=upf)

    uc0 = Constant(vel)
    u_conv = Constant([vel])

    n = FacetNormal(mesh)
    flux_nU = u*(dot(u_conv, n) + abs(dot(u_conv, n)))/2
    flux = flux_nU('+') - flux_nU('-')
    bflux = flux_nU('+')
    eq = Constant(1/dt)*(u - upf)*v*dx - u*uc0*v.dx(0)*dx + flux*jump(v)*dS + bflux*v*ds

    a, L = lhs(eq), rhs(eq)


    As.append(assemble(a).array())
    
print As[0]
print As[1]

print numpy.linalg.norm(As[0]-As[1])
