from dolfin import *

mesh = UnitSquareMesh(10, 10)

V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

x = Expression('x[0]', degree=1)
y = Expression('x[1]', degree=1)

eq = inner(grad(u)+as_vector((x, y)), grad(v))*dx
a, L = system(eq)
print assemble(a)
print assemble(L)

a = inner(grad(u) + u*as_vector((x, y)), grad(v))*dx
print assemble(a)
