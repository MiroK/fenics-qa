from dolfin import *
''' Eigenvalue problem
    -laplace(u) = lmbda * u on (0, 1)
    u(0) = u(1) = 0
'''


mesh = UnitIntervalMesh(100)
V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

# Lhs mult. by v and integrated by-parts gives
a = inner(grad(u), grad(v))*dx
# Rhs mult by v and integrated gives
m = inner(u, v)*dx
# Auxiliary
L = Constant(0)*v*dx

bc = DirichletBC(V, Constant(0.), DomainBoundary())

# Assemble matrices
A = PETScMatrix()
M = PETScMatrix()
b = PETScVector()
assemble_system(a, L, bc, A_tensor=A, b_tensor=b)
assemble_system(m, L, bc, A_tensor=M, b_tensor=b)

# Setup the eigensolver
eigensolver = SLEPcEigenSolver(A, M)
eigensolver.solve()

# Get the eigenpairs real/complex eigenvalues, real/complex eigenvectors
for i in range(eigensolver.get_number_converged()):
    r, c, rx, cx = eigensolver.get_eigenpair(i)
    print r
