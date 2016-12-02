from dolfin import *
import sys

n = int(sys.argv[1])
mesh = UnitSquareMesh(n, n)
# V = VectorFunctionSpace(mesh, 'CG', 1)
V = FunctionSpace(mesh, 'N1curl', 2)

f = Expression(('x[1]', '-x[0]'), degree=1)
f_e0 = Expression(('x[1]', '0'), degree=1)

u = interpolate(f, V)
e0 = interpolate(Constant((1, 0)), V)

u_e0 = (dot(u, e0)/dot(e0, e0))*e0
u_e0 = project(u_e0, V)

print errornorm(f_e0, u_e0, 'L2')
