"""
FEniCS tutorial demo program:
Nonlinear Poisson equation with Dirichlet conditions
in x-direction and homogeneous Neumann (symmetry) conditions
in all other directions. The domain is the unit hypercube in
of a given dimension.

-div(q(u)*grad(u)) = 0,
u = 0 at x=0, u=1 at x=1, du/dn=0 at all other boundaries.
q(u) = (1+u)^m

Solution method: automatic, i.e., by a NonlinearVariationalProblem/Solver
(Newton method), with automatic UFL computation of the derivative.
"""
from dolfin import *
import numpy, sys

# Usage:   ./vp2_np.py m|a |g|l degree nx ny nz
# Example: ./vp2_np.py m    l   1      3  4
J_comp = 'a'  # m (manual) or a (automatic) computation of J
answer = 'l'  # g (GMRES) or l (sparse LU) solver
iterative_solver = True if answer == 'g' else False

# Create mesh and define function space
degree = 1
divisions = [8]
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
v  = TestFunction(V)
du = TrialFunction(V)
u_ = Function(V)  # most recently computed solution
F  = inner(q(u_)*grad(u_), grad(v))*dx

# J must be a Jacobian (Gateaux derivative in direction of du)
if J_comp == 'm':
    J = inner(q(u_)*grad(du), grad(v))*dx + \
        inner(Dq(u_)*du*grad(u_), grad(v))*dx
else:
    J = derivative(F, u_, du)

# Compute solution
problem = NonlinearVariationalProblem(F, u_, bcs, J)
