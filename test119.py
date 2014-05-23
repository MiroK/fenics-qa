from dolfin import *

mesh = RectangleMesh(0.0, 0.0, 1.0, 1.0, 100, 100)
V = FunctionSpace(mesh, "Lagrange", 1)
V_vec = VectorFunctionSpace(mesh, "Lagrange", 1)

c = project(Expression('x[0]'), V)

t = Timer('foo')
t.start()
for i in range(0):
    v = as_vector((1, 2))
    d = project(c*v,V_vec)
t.stop()

t = Timer('bar')
t.start()
for i in range(1):
    d = Function(V_vec)
    x = interpolate(Constant((1, 2)), V_vec)
    c.vector().inner(x.vector())
t.stop()

print 'foo', timing('foo')


