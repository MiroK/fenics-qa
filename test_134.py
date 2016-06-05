from dolfin import *

def Max(a, b): return (a+b+abs(a-b))/Constant(2)
def Min(a, b): return (a+b-abs(a-b))/Constant(2)

mesh = RectangleMesh(Point(-1, -1), Point(1, 1), 10, 10, 'crossed')
V = VectorFunctionSpace(mesh, 'CG', 1)
v = interpolate(Expression(('x[0]', '-x[1]'), degree=1), V)

tval = sqrt(dot(v, v))
fval = Constant(0.5)

P = FunctionSpace(mesh, 'DG', 0)
pmax = project(Max(tval, fval), P)
pmin = project(Min(tval, fval), P)

plot(pmax, title='max')
plot(pmin, title='min')
interactive()
