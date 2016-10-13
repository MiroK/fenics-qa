from dolfin import *

# Create mesh
mesh = UnitSquareMesh(32, 32)

# Define function spaces and mixed (product) space
TCG = TensorFunctionSpace(mesh, "CG", 2)
VCG = VectorFunctionSpace(mesh, "CG", 1)
W = TCG * VCG

nu = 0.10
E = 1.0
mu = Constant(E / (2.0 * (1.0 + nu)))
lmbda = Constant(E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)))

def compliance(sigma):
    return 1.0 / 2.0 * mu * (sigma - lmbda / (2 * mu + 2 * lmbda) * tr(sigma) * Identity(2))

# Define trial and test functions
(sigma, u) = TrialFunctions(W)
(tau, v) = TestFunctions(W)

# Define source function
f = Constant((0.0, 0.0))

# Define variational form
a = (inner(compliance(sigma), tau) + inner(div(tau), u) + inner(div(sigma), v))*dx
L = - inner(f, v)*dx# + inner(u0, dot(n, tau))*ds(1) + inner(u0, dot(n, tau))*ds(2)

# Define function G such that G \cdot n = g
class BoundarySource(Expression):
    def __init__(self, mesh):
        self.mesh = mesh
    def eval_cell(self, values, x, ufc_cell):
        cell = Cell(self.mesh, ufc_cell.index)
        n = cell.normal(ufc_cell.local_facet)
        values[0] = 1.0
        values[1] = 0.0
        values[2] = 0.0
        values[3] = 1.0
    def value_shape(self):
        return (2, 2)

G = BoundarySource(mesh)

# Define essential boundary
def boundary(x):
    return x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS

bc = DirichletBC(W.sub(0), G, boundary)

# Compute solution
w = Function(W)
solve(a == L, w, bc)
(sigma, u) = w.split()

plot(u)
interactive()
