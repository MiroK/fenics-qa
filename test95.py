from dolfin import *

mesh = UnitSquareMesh(25, 25)
V = FunctionSpace(mesh, 'CG', 1)

foo = Expression('sin(2*pi*x[0])*sin(2*pi*x[1])*pow(t, 2)', t=0.3)
u0 = interpolate(foo, V)

foo.t=0.6
u = interpolate(foo, V)

# v = u1*(1-alpha) + u2*alpha
v = Function(V)
alpha = 0.5
v.vector()[:] = (1 - alpha)*u0.vector().array() + alpha*u.vector().array()

plot(v, interactive=True)

