from dolfin import *

f = Expression('(2*pi*pi+1)*sin(pi*x[0])*sin(pi*x[1])+1.0')
g = Constant('1.0')

mesh = UnitSquareMesh(20, 20)

V = FunctionSpace(mesh, 'CG', 1)
M = MixedFunctionSpace([V, V])
u, ambda = TrialFunctions(M)
v, mu = TestFunctions(M)

a = inner(nabla_grad(u), nabla_grad(v))*dx + u*v*dx + (ambda*v + mu*u)*ds
L = f*v*dx + mu*g*ds

A, b = assemble_system(a, L)
A.ident_zeros()

u_ambda = Function(M)
solve(A, u_ambda.vector(), b)
u, ambda = u_ambda.split()

plot(u, interactive=True)
plot(ambda, interactive=True)
