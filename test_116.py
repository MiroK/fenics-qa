from dolfin import *

# Elasticity parameters
E, nu = 10.0, 0.3
mu, lambda_ = Constant(E/(2*(1 + nu))), Constant(E*nu/((1 + nu)*(1 - 2*nu)))# Strain-rate

# Strain
epsilon = lambda u: sym(grad(u))

# Stress
sigma = lambda u: 2*mu*epsilon(u) + lambda_*tr(epsilon(u))*Identity(2)

# Nullspace
nullspace = [Constant((1, 0)), Constant((0, 1)), Expression(('-x[1]', 'x[0]'), degree=1)]


mesh = UnitSquareMesh(100, 100)
V = VectorFunctionSpace(mesh, 'CG', 1)       # Space for displacement


# Let's first show that the energy form without constraints is not pos.def
u, v = TrialFunction(V), TestFunction(V)
# Energy
a = inner(sigma(u), epsilon(v))*dx
A = PETScMatrix()
assemble(a, A)
# Every z is eigenv with eigenvalue 0
for z in nullspace:
    v = interpolate(z, V).vector()
    null = v.copy()
    print null.norm('l2'),
    A.mult(v, null)
    print null.norm('l2')
print
##############################################################################

# Confirm that A has three zero eigenvalue with slepc
eigensolver = SLEPcEigenSolver(A)
eigensolver.parameters['problem_type'] = 'hermitian'
eigensolver.parameters['spectrum'] = 'smallest magnitude'
eigensolver.solve(5)

assert eigensolver.get_number_converged() > 3
for i in range(4):
    w = eigensolver.get_eigenvalue(i)
    print w
print
##############################################################################

# Get rid of singularity
M = VectorFunctionSpace(mesh, 'R', 0, 3)     # Space for all multipliers
W = MixedFunctionSpace([V, M])
u, mus = TrialFunctions(W)
v, nus = TestFunctions(W)

# Energy
a = inner(sigma(u), epsilon(v))*dx

# Lagrange multipliers contrib to a
for i, e in enumerate(nullspace):
    mu = mus[i]
    nu = nus[i]
    a += mu*inner(v, e)*dx + nu*inner(u, e)*dx

# See that there are now no zero eigenvalues
A = PETScMatrix()
assemble(a, A)

eigensolver = SLEPcEigenSolver(A)
eigensolver.parameters['problem_type'] = 'hermitian'
eigensolver.parameters['spectrum'] = 'smallest magnitude'
eigensolver.solve(5)

assert eigensolver.get_number_converged() > 1
for i in range(1):
    w = eigensolver.get_eigenvalue(i)
    print w
##############################################################################
