from dolfin import *

a=5.e-2; #cube size
# Sub domain for Periodic boundary condition
class PeriodicBoundary(SubDomain):

    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        # return True if on left or bottom boundary AND NOT on one of the two slave edges
        return bool(( near(x[0], -a/2.) or near(x[1], -a/2.) or near(x[2], -a/2.)) and 
            (not ((near(x[0], a/2.) and near(x[2], -a/2.)) or 
                  (near(x[0], -a/2.) and near(x[2], a/2.)) or
                  (near(x[1], a/2.) and near(x[2], -a/2.))or
                  (near(x[1], -a/2.) and near(x[2], a/2.)))) and on_boundary)

    # Map right boundary (H) to left boundary (G)
    def map(self, x, y):

        if near(x[0], a/2.) and near(x[2], a/2.):
            y[0] = x[0] - a/2
            y[1] = x[1] 
            y[2] = x[2] - a/2
        elif near(x[1], a/2.) and near(x[2], a/2.):
            y[0] = x[0] 
            y[1] = x[1] - a/2
            y[2] = x[2] - a/2
        elif near(x[0], a/2):
            y[0] = x[0] - a/2
            y[1] = x[1]
            y[2] = x[2]
        elif near(x[1], a/2):
            y[0] = x[0]
            y[1] = x[1] - a/2
            y[2] = x[2] 
        elif near(x[2], a/2):
            y[0] = x[0]
            y[1] = x[1]
            y[2] = x[2] - a/2



mesh = BoxMesh(Point(-a/2, -a/2, -a/2), Point(a/2, a/2, a/2), 36,36,36)

V = VectorElement("Lagrange", mesh.ufl_cell(), 1)
W = MixedElement([V, V])
Vcomplex = FunctionSpace(mesh, W, constrained_domain=PeriodicBoundary())

print Vcomplex.dim(), FunctionSpace(mesh, W).dim()
