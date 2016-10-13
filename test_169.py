from dolfin import *
import numpy as np
import sys

# Right end point 
xmax = 1.
ymax = 1.

# Left end point 
xmin = -1.
ymin = -1.

# Number of elements
nel_x = 100
nel_y = 100

# Define mesh
mesh = RectangleMesh(Point(xmin, ymin), Point(xmax, ymax), \
                     nel_x, nel_y, "crossed")
#plot(mesh, interactive=True)

# Define spatial tolerance
tol_x = (xmax - xmin) / nel_x / 2.
tol_y = (ymax - ymin) / nel_y / 2.

# Create classes for defining parts of the boundaries and the interior
# of the domain
class Left(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], xmin, tol_x) and on_boundary

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], xmax, tol_x) and on_boundary

# Elasticity parameters
nu = float(sys.argv[1])
E = 1.0
mu = Constant(E / (2.0 * (1.0 + nu)))
lmbda = Constant(E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)))

# Strain
def eps(v):
    return sym(grad(v))

# Stress
def sigma(v):
    return 2.0 * mu * eps(v) + lmbda * tr(eps(v)) * Identity(len(v))

# Define Dirichlet boundary (0 at (x,y) = (0,0))
def boundary(x):
    return (near(x[0], 0.0, tol_x) and near(x[1], 0.0, tol_y))  # for reasons
           # of stability, otherwise the rectangle is rotated

# Initialize sub-domain instances
left = Left()
right = Right()

# Initialize mesh function for boundary domains
boundaries = FacetFunction("size_t", mesh)
boundaries.set_all(0)
left.mark(boundaries, 1)
right.mark(boundaries, 2)
ds = Measure("ds")[boundaries]

# Define active stress at boundaries
g_L = Constant((-0.001, 0.0))
g_R = Constant((0.001, 0.0))

# Define function space and basis functions
V = VectorFunctionSpace(mesh, "CG", 2)
u = TrialFunction(V)
v = TestFunction(V)

# Define boundary condition
u0 = Constant((0.0, 0.0))
bc = DirichletBC(V, u0, boundary, method="pointwise")

# Define new measures associated with the interior domains and
# exterior boundaries
ds = Measure("ds")[boundaries]

# Define variational form
a = inner(sigma(u), sym(grad(v))) * dx
L = inner(g_L, v) * ds(1) + inner(g_R, v) * ds(2)

# Solve problem
u = Function(V)
problem = dolfin.LinearVariationalProblem(a, L, u, bc)
solver = dolfin.LinearVariationalSolver(problem)
solver.parameters["linear_solver"] = "default"
solver.parameters["preconditioner"] = "ilu"
solver.solve()

## Save nodal positions
#fid = File("results/solution.pvd")
#u.rename("u", "u")
#fid << u

# Plot solution and gradient
#plot(u, interactive=True)
print "Theoretical Poisson's ratio: ", nu
print "Simulated Poisson's ratio: ", \
    - np.log(1.0 - abs(u(0., 1.)[1])) / np.log(1.0 + abs(u(1., 0.)[0]))
