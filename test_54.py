from dolfin import *

parameters['form_compiler']['representation'] = 'quadrature'
f = Expression('std::pow(std::sqrt(x[0]*x[0]+x[1]*x[1]), -1./6)', degree=4)

mesh = RectangleMesh(-1, -1, 1, 1, 101, 101)
V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

a = inner(grad(u), grad(v))*dx + inner(u, v)*dx

#W = FunctionSpace(mesh, 'DG', 0)
#f = project(f, W)
L = inner(f, v)*dx


u = Function(V)
solve(a == L, u)

plot(u)
interactive()
