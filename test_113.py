from dolfin import *

mesh = UnitSquareMesh(20, 20)
V = FunctionSpace(mesh, 'CG', 2)
f = Expression("x[0]*x[0] + x[1]", degree=2)

v = TestFunction(V)
b = assemble(inner(f, v)*dx)
fv = project(f, V, solver='lu')

b2 = assemble(inner(fv, v)*dx)

b2.axpy(-1, b)
print b2.norm('linf')

fv2 = project(fv, V, solver='lu').vector()
fv2.axpy(-1, fv.vector())
print fv2.norm('linf')
