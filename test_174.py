from dolfin import *
import numpy as np

#        /  f3 if f1 > f2
# f_4 =  |
#        \  f1 otherwise

mesh = IntervalMesh(1000, -1, 1)
V = FunctionSpace(mesh, 'CG', 1)
f1 = interpolate(Expression('1-x[0]*x[0]', degree=2), V)
f2 = interpolate(Constant(3/4.), V)
f3 = interpolate(Expression('3*x[0]*x[0]', degree=2), V)

# By Projection
f4 = project(conditional(gt(f1, f2), f3, f1), V)

# By manipulating vectors
f5 = Function(V, f1.vector())

F5 = f5.vector().get_local()
F2 = f2.vector().get_local()
F3 = f3.vector().get_local()

gt = np.where(F5 > F2)[0]
F5[gt] = F3[gt]

f5.vector().set_local(F5)
f5.vector().apply('insert')

print '|f1-f2|_2 =', sqrt(assemble(inner(f1-f5, f1-f5)*dx))

plot(f4)
interactive()
