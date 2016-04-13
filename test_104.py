from dolfin import *

parameters["linear_algebra_backend"] = "PETSc"

mesh = BoxMesh(Point(-1, -1, -1), Point(1, 1, 1), 10, 10, 10)

# Elasticity parameters
E = 1.0e9
nu = 0.3
mu = E/(2.0*(1.0 + nu))
lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu))

# Load
omega = 300.0
rho = 10.0
f = Expression(("rho*omega*omega*x[0]", "rho*omega*omega*x[1]", "0.0"),
               omega=omega, rho=rho, degree=2)

# Variational formulation
V = VectorFunctionSpace(mesh, "Lagrange", 1)
# Stress computation
def sigma(v):
    return 2.0*mu*sym(grad(v)) + lmbda*tr(sym(grad(v)))*Identity(len(v))

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
a = inner(sigma(u), grad(v))*dx
L = inner(f, v)*dx

# Set up boundary condition on inner surface
c = Constant((0.0, 0.0, 0.0))
bc = DirichletBC(V, c, 'near(x[0], -1)')

A = 8*sqrt(10)/3
# Nullspace of righd motions
Z = [Constant((1/8., 0, 0)),
     Constant((0, 1/8., 0)),
     Constant((0, 0, 1/8.)),
     Expression(('0', 'x[2]/A', '-x[1]/A'), A=A),
     Expression(('-x[2]/A', '0', 'x[0]/A'), A=A),
     Expression(('x[1]/A', '-x[0]/A', '0'), A=A)]
# As vectors
Z = [interpolate(z, V).vector() for z in Z]
# For FEniCS
for zi in Z:
    for zj in Z:
        print zi.inner(zj),
    print
Z = VectorSpaceBasis(Z)

# Assemble system, applying boundary conditions and preserving
# symmetry)
A, b = assemble_system(a, L, bc)

# Create solution function
u = Function(V)

# Create PETSC smoothed aggregation AMG preconditioner and attach near
# null space
pc = PETScPreconditioner("petsc_amg")
pc.set_nullspace(Z)
Z.orthogonalize(b)

# Create CG Krylov solver and turn convergence monitoring on
solver = PETScKrylovSolver("cg", pc)
solver.parameters["monitor_convergence"] = True
solver.parameters["relative_tolerance"] = 1E-8

# Set matrix operator
solver.set_operator(A);

# Compute solution
solver.solve(u.vector(), b);

# Plot solution
plot(u, interactive=True)
