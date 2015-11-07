from sympy import symbols, lambdify, besselj, I
from scipy.special import jv as scipy_besselj
from dolfin import *

mesh = UnitCubeMesh(20, 20, 20)

# Space where the exact solution is represented
Ve = FunctionSpace(mesh, 'DG', 4)
dof_x = Ve.dofmap().tabulate_all_coordinates(mesh).reshape((-1, 3))
X, Y, Z = dof_x[:, 0], dof_x[:, 1], dof_x[:, 2]

# Suppose the solution is besselj(1, x[0]) + besselj(2, x[1]) + besselj(3, x[2])
x, y, z = symbols('x, y, z')
u = (0.00073 - 0.00052*I)*besselj(0,(1.84 - 1.84*I)*x)
u_lambda = lambdify(x, u, ['numpy', {'besselj': scipy_besselj}])
# Interpolate manually
# timer = Timer('eval')
# timer.start()
# u_values = u_lambda(X, Y, Z)
# print 'Evaluated %d dofs in % s' % (Ve.dim(), timer.stop())
# Exact solution in Ve
# u = Function(Ve)
# u.vector()[:] = u_values

# Solution space and the hypot. numerical solution
# V = FunctionSpace(mesh, 'CG', 1)
# uh = interpolate(Constant(1), V)

# print errornorm(u, uh)

print u_lambda(1.0)
