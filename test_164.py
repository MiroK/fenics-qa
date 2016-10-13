from dolfin import *

mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, 'CG', 1)
v = TestFunction(V)
f0 = Function(V)
f1 = Function(V)

L0 = inner(f0, v)*dx
L1 = inner(f1, v)*dx

print L0.equals(L1)
L0 = replace(L0, {f0: f1})
print L0.equals(L1)
