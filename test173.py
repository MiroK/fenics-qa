from dolfin import *

mesh = UnitSquareMesh(10, 10)

V = VectorFunctionSpace(mesh, 'CG', 2)
Q = FunctionSpace(mesh, 'CG', 1)
W = MixedFunctionSpace([V, Q])

f = Expression(('0', '1', 'x[0]+x[1]'))
up = interpolate(f, W)
u, p = up.split()

plot(u, title='u')
plot(p, title='p')
interactive()
