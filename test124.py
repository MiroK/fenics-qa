from dolfin import *

mesh = UnitIntervalMesh(10)
V = FunctionSpace(mesh, 'CG', 1)

du = TrialFunction(V)
u = Function(V)
v = TestFunction(V)

F = inner(u**3, v)*dx
dF = derivative(F, u, u)

print assemble(dF).array()

uk = interpolate(Expression('x[0]'), V)
dF_uk = action(dF, uk)

foo = assemble(dF_uk)
print foo.array()
bar = Function(V, foo)

baz = interpolate(Expression('3*x[0]*x[0]'), V)
qux = assemble(inner(baz, v)*dx)
print qux.array()



