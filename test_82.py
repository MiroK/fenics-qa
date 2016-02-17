from dolfin import *

mesh = UnitSquareMesh(3, 3)
V = FunctionSpace(mesh, 'CG', 1)
W = MixedFunctionSpace([V, V])

v1 = interpolate(Expression('x[0]'), V)
print v1.vector().norm('l2')

v2 = interpolate(Expression('2*x[0]'), V)
print v2.vector().norm('l2')
v1.assign(v2)
print v1.vector().norm('l2')

w = Function(W)
assign(w.sub(0, deepcopy=False), v1)
print w.vector().norm('l2')
