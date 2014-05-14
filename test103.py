from dolfin import *

# Define mesh and function space
mesh = RectangleMesh(-10, -10, 10, 10, 20, 20)
V = FunctionSpace(mesh, 'Lagrange', 1)

# Apply essential boundary conditions
bc = DirichletBC(V, Constant(0.), DomainBoundary())

# Define functions
u = TrialFunction(V)
v = TestFunction(V)

Pot = Expression('fabs(x[0]+x[1])+fabs(x[0]-x[1])')

# Define problem
a = inner(grad(u), grad(v))*dx + Pot*u*v*dx
m = u*v*dx
L = Constant(0.)*v*dx

# Assemble stiffness matrix
A = PETScMatrix()
M = PETScMatrix()
_ = PETScVector()

assemble_system(a, L, bc, A_tensor=A, b_tensor=_)
assemble_system(m, L, bc, A_tensor=M, b_tensor=_)

# Create eigensolver
eigensolver = SLEPcEigenSolver(A, M)
eigensolver.parameters['spectrum'] = 'smallest magnitude'
eigensolver.parameters['solver']   = 'lapack'
eigensolver.parameters['tolerance'] = 1.e-15

# Solve for eigenvalues
eigensolver.solve()

u = Function(V)
for i in range(0, eigensolver.get_number_converged()):
    # Extract next eigenpair
    r, c, rx, cx = eigensolver.get_eigenpair(i)
    if r > 1:
        print 'eigenvalue:', r

        # Assign eigenvector to function
        u.vector()[:] = rx

        # Plot eigenfunction
        plot(u)
        interactive()
