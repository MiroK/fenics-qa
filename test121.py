from dolfin import *

J_comp = 'm'
answer = 'l'
iterative_solver = True if answer == 'g' else False

# Create mesh and define function space
degree = 1
divisions = 10, 5
d = len(divisions)
domain_type = [UnitIntervalMesh, UnitSquareMesh, UnitCubeMesh]
mesh = domain_type[d-1](*divisions)
V = FunctionSpace(mesh, 'Lagrange', degree)


# Define boundary conditions
tol = 1E-14
def left_boundary(x, on_boundary):
    return on_boundary and abs(x[0]) < tol

def right_boundary(x, on_boundary):
    return on_boundary and abs(x[0]-1) < tol

Gamma_0 = DirichletBC(V, Constant(0.0), left_boundary)
Gamma_1 = DirichletBC(V, Constant(1.0), right_boundary)
bcs = [Gamma_0, Gamma_1]

# Choice of nonlinear coefficient
m = 2

def q(u):
    return (1+u)**m

def Dq(u):
    return m*(1+u)**(m-1)

# Define variational problem
v = TestFunction(V)
u = TrialFunction(V)
F = inner(q(u)*nabla_grad(u), nabla_grad(v))*dx
u_ = interpolate(Constant(10), V) # most recently computed solution
F = action(F, u_)

# J must be a Jacobian (Gateaux derivative in direction of du)
if J_comp == 'm':
    J = inner(q(u_)*nabla_grad(u), nabla_grad(v))*dx + \
        inner(Dq(u_)*u*nabla_grad(u_), nabla_grad(v))*dx
else:
    J = derivative(F, u_, u)

# Compute solution

u_ = interpolate(Constant(10), V) # most recently computed solution
problem = NonlinearVariationalProblem(F, u_, bcs, J)
solver = NonlinearVariationalSolver(problem)
prm = solver.parameters
info(prm, True)
prm['newton_solver']['absolute_tolerance'] = 1E-8
prm['newton_solver']['relative_tolerance'] = 1E-7
prm['newton_solver']['maximum_iterations'] = 25
prm['newton_solver']['relaxation_parameter'] = 1.0
if iterative_solver:
    prec = 'jacobi' if 'jacobi' in zip(*krylov_solver_preconditioners())[0] else 'ilu'
    prm['newton_solver']['linear_solver'] = 'gmres'
    prm['newton_solver']['preconditioner'] = prec
    prm['newton_solver']['krylov_solver']['absolute_tolerance'] = 1E-9
    prm['newton_solver']['krylov_solver']['relative_tolerance'] = 1E-7
    prm['newton_solver']['krylov_solver']['maximum_iterations'] = 1000
    prm['newton_solver']['krylov_solver']['monitor_convergence'] = True
    prm['newton_solver']['krylov_solver']['nonzero_initial_guess'] = False
    prm['newton_solver']['krylov_solver']['gmres']['restart'] = 40
    prm['newton_solver']['krylov_solver']['preconditioner']['structure'] = 'same_nonzero_pattern'
    prm['newton_solver']['krylov_solver']['preconditioner']['ilu']['fill_level'] = 0
    prm['newton_solver']['krylov_solver']['nonzero_initial_guess'] = True
PROGRESS = 16
set_log_level(PROGRESS)
solver.solve()

print """
Solution of the nonlinear Poisson problem div(q(u)*nabla_grad(u)) = f,
with f=0, q(u) = (1+u)^m, u=0 at x=0 and u=1 at x=1.
%s
""" % mesh

# Find max error
u_exact = Expression('pow((pow(2, m+1)-1)*x[0] + 1, 1.0/(m+1)) - 1', m=m)
u_e = interpolate(u_exact, V)
import numpy
diff = numpy.abs(u_e.vector().array() - u_.vector().array()).max()
print 'Max error:', diff

plot(u_e, interactive=True)
