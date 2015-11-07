from dolfin import *

f = Expression('cos(2*pi*x[0])*x[1]')

mesh = UnitSquareMesh(3, 3)

# Sub domain for Periodic boundary condition
class PeriodicBoundary(SubDomain):
    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        return bool(x[0] < DOLFIN_EPS and x[0] > -DOLFIN_EPS and on_boundary)

    # Map right boundary (H) to left boundary (G)
    def map(self, x, y):
        y[0] = x[0] - 1.0
        y[1] = x[1]

# Constrained space
Vp = FunctionSpace(mesh, 'CG', 2, constrained_domain=PeriodicBoundary())
p = interpolate(f, Vp)
print p.vector().array()

# Unconstrained space
V = FunctionSpace(mesh, 'CG', 2)
v = interpolate(p, V)
print v.vector().array()
