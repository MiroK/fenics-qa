from dolfin import *

# define periodic bc
class PeriodicBC(SubDomain):
      def __init__(self, tolerance=DOLFIN_EPS, length = 1.0):
          SubDomain.__init__(self)
          self.tol = tolerance
          self.length = length
      def inside(self, x, on_boundary):
          return bool((near(x[0], 0) or near(x[1], 0) and on_boundary) and \
                      (not ((near(x[0], 0) and near(x[1], self.length)) or \
                      (near(x[0], self.length) and near(x[1], 0)))) and on_boundary)
      def map(self, x, y):
          L = self.length
          if near(x[0], L) and near(x[1], L):
             y[0] = x[0] - L
             y[1] = x[1] - L
          elif near(x[0], L):
             y[0] = x[0] - L
             y[1] = x[1]
          else:  
             y[0] = x[0]
             y[1] = x[1] - L


# create mesh
mesh = RectangleMesh(Point(0.0,0.0),Point(1.0,1.0),4,4,"left/right")

# build function space
PBC = PeriodicBC(length = 1.0)
V  = FunctionSpace(mesh, "CG", 1, constrained_domain=PBC)
W  = FunctionSpace(mesh, "CG", 1, constrained_domain=PBC)          
MS = MixedFunctionSpace([V,W])
print MS.dim()

V = FiniteElement('Lagrange', mesh.ufl_cell(), 1)
W = MixedElement([V, V])
W = FunctionSpace(mesh, W, constrained_domain=PBC)          
print W.dim()
