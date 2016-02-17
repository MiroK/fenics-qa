from dolfin import *

# Setup problem
mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, 'CG', 1)
bc = DirichletBC(V, 1.0, 'on_boundary')

u = Function(V)
v = TestFunction(V)

x = SpatialCoordinate(mesh)
r = sqrt(x[0]**2 + x[1]**2)
rho = 1.0/r**3
F = 8*inner(grad(u), grad(v))*dx + rho * inner(u**5, v)*dx  +\
    (-1.0/8.0)*inner(u, v)*dx

du = TrialFunction(V)
J = derivative(F, u, du)

params = {'nonlinear_solver': 'snes',
          'snes_solver': {'linear_solver': 'lu',
          'maximum_iterations': 100,
          'sign': 'nonnegative',
          'report': True}}

# Setup solve
u.interpolate(Constant(-1000.0))
problem = NonlinearVariationalProblem(F, u, bc, J=J)
solver  = NonlinearVariationalSolver(problem)
solver.parameters.update(params)
info(solver.parameters, True)

# Solve the problem
(iter, converged) = solver.solve()

# Get SNES - and this does not work :(
snes = solver.snes_solver().snes()
print type(snes)
