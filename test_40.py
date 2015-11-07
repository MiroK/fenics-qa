from dolfin import *

mesh = UnitSquareMesh(20, 20)

V = FunctionSpace(mesh, 'CG', 1)
v = Function(V).vector()
print type(v)

v = as_backend_type(v)
print type(v)
print 'length of v', v.vec().size
