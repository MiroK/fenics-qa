from dolfin import *
mesh = UnitSquareMesh(10,10)
TS = TensorFunctionSpace(mesh, 'CG', 1)
Q = interpolate(Constant(((1, 2), (3, 4))), TS)

u = project(dot(Q,Q), TS) # works
v = project(as_tensor(Q[i, k]*Q[k, j], (i, j)), TS)

u.vector().axpy(-1, v.vector())
print u.vector().norm('linf')
print v.vector().norm('linf')
