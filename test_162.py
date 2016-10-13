from dolfin import *

import math
import ufl

import numpy as np

# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True
ffc_options = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
               "precompute_ip_const": True}

# User inputs
nel = 100

q = Constant(1.0);


# Create mesh and define function space
mesh = IntervalMesh(nel,0.0,1.0)

V = FunctionSpace(mesh, "Lagrange", 1) # Linear Lagrange Fintie Elements on 1D Mesh

# Mark Boundaries
left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)

# Define Dirichlet boundary (u = 0 at x = 0)

bcl = DirichletBC(V, Constant(0.0), left)

# Define functions
du = TrialFunction(V)            # Incremental Displacement
v  = TestFunction(V)             # Test function
u  = Function(V)                 # Displacement from previous iteration

# Loading Looping


# Total Potential Energy
W =  (1.0/3.0) * (u ** 3) * dx(mesh) - q * u * dx(mesh)

# Compute first variation of Pi (directional derivative about u in the direction of v)
F = derivative(W, u, v)

# Compute Jacobian of F
J = derivative(F, u, du)

# Solve variational problem
problem = NonlinearVariationalProblem(F, u, bcl, J,
                                      form_compiler_parameters=ffc_options)
solver  = NonlinearVariationalSolver(problem)

prm = solver.parameters
prm['newton_solver']['absolute_tolerance'] = 1E-8
prm['newton_solver']['relative_tolerance'] = 1E-7
prm['newton_solver']['maximum_iterations'] = 100
prm['newton_solver']['relaxation_parameter'] = 2.

solver.solve()
