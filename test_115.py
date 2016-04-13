from mshr import *
from dolfin import *

#set_log_level(1)

# Test for PETSc
if not has_linear_algebra_backend("PETSc"):
    print("DOLFIN has not been configured with PETSc. Exiting.")
    exit()
    
# Do it parallel
parameters["num_threads"] = 4;

# Set backend to PETSC
parameters["linear_algebra_backend"] = "PETSc"

class FreeBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

class LeftBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and abs(x[0] + 2) < DOLFIN_EPS

class RightBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and abs(x[0] - 2) < DOLFIN_EPS


def build_nullspace(V, x):
    #Function to build null space for 2D elasticity
    #
    # Create list of vectors for null space
    nullspace_basis = [x.copy() for i in range(2)]

    # Build translational null space basis - 2D problem
    #
    V.sub(0).dofmap().set(nullspace_basis[0], 1.0);
    V.sub(1).dofmap().set(nullspace_basis[1], 1.0);

    for x in nullspace_basis:
        x.apply("insert")

    r = interpolate(Expression(('-x[1]', 'x[0]')), V).vector()
    nullspace_basis.append(r)

    # Create vector space basis and orthogonalize
    basis = VectorSpaceBasis(nullspace_basis)
    basis.orthonormalize()

    return basis

# Initialize subdomain instances
free_b = FreeBoundary()
left_b = LeftBoundary()
right_b = RightBoundary()
    
# Create geo and mesh
res = 100
r1 = Rectangle(Point(-2., -2.), Point(2., 2.))
r2 = Rectangle(Point(-0.5, -0.5), Point(0.5, 0.5))
geo = r1-r2
mesh = generate_mesh(geo,res)

# Initialize mesh function for boundary domains
boundaries = FacetFunction("size_t", mesh)
boundaries.set_all(0)
free_b.mark(boundaries, 1)
left_b.mark(boundaries, 2)
right_b.mark(boundaries, 3)

# Create function space
V = VectorFunctionSpace(mesh, 'Lagrange', 1)

# Create test and trial functions, and source term
u = TrialFunction(V)
v = TestFunction(V)

# Define new measures associated with the exterior boundaries
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

# Elasticity parameters
E, nu = 10, 0.3
mu, lmbda = E/(2.0*(1.0 + nu)), E*nu/((1.0 + nu)*(1.0 -2.0*nu))

# Stress relation
I = Identity(2)
epsilon = sym(nabla_grad(u))
sigma = 2*mu*epsilon + lmbda*tr(epsilon)*I

# Boundary conditions - Neumann (stress on L&R side)
gL = Constant((-1.0, 0.0))
gR = Constant((1.0, 0.0))

# Governing balance equation
a = inner(sigma, sym(nabla_grad(v)))*dx
L = dot(gL, v)*ds(2) + dot(gR, v)*ds(3)

A, b = assemble_system(a, L)
print "Norm(b)=",b.norm("l2") # Check 

# Set up PDE and solve
u = Function(V)

# AMG). The solution vector is passed so that it can be copied to
# generate compatible vectors for the nullspace.
null_space = build_nullspace(V, u.vector())

# Attach near nullspace to matrix and orthogonalize b
as_backend_type(A).set_near_nullspace(null_space)
null_space.orthogonalize(b) 

# Create PETSC smoothed aggregation AMG preconditioner and attach near
# null space
pc = PETScPreconditioner("petsc_amg")

# Use Chebyshev smoothing for multigrid
PETScOptions.set("mg_levels_ksp_type", "chebyshev")
PETScOptions.set("mg_levels_pc_type", "jacobi")

# Improve estimate of eigenvalues for Chebyshev smoothing
PETScOptions.set("mg_levels_esteig_ksp_type", "cg")
PETScOptions.set("mg_levels_ksp_chebyshev_esteig_steps", 50)

# Create CG Krylov solver and turn convergence monitoring on
solver = PETScKrylovSolver("cg", pc)
solver.parameters["monitor_convergence"] = True

# Set matrix operator
solver.set_operator(A);

# Compute solution
solver.solve(u.vector(), b);

plot(u, mode = "displacement", interactive=True, scalarbar=True, title="u")
