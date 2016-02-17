# Begin demo
from dolfin import *

# Create mesh and define function space
mesh = UnitSquareMesh(32, 32)
V = FunctionSpace(mesh, "Lagrange", 1)

# Exterior domain (PUT ATTENTION HERE!)
class Exterior(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > 0.75 - DOLFIN_EPS or x[0] < 0.25 + DOLFIN_EPS or x[1] > 0.75 - DOLFIN_EPS or x[1] < 0.25 + DOLFIN_EPS

# Define Dirichlet boundary (x = 0 or x = 1)
def boundary(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

# Define boundary condition (AND HERE!)
exterior = Exterior()
u0 = Constant(0.0)
bc1 = DirichletBC(V, u0, boundary)
bc2 = DirichletBC(V, u0, exterior)
bcs = [bc1, bc2]

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)")
g = Expression("sin(5*x[0])")
a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*ds

# Compute solution
u = Function(V)
solve(a == L, u, bcs)

# Plot solution
plot(u, interactive=True)
