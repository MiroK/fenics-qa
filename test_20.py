from dolfin import *
import numpy as np

fe = Expression('k*sin(k*x[0])', k=0)

mesh = UnitIntervalMesh(100)
V = FunctionSpace(mesh, 'CG', 1)
f = interpolate(fe, V)

for k in np.linspace(1, 2, 100):
    fe.k = k
    f.vector()[:] = interpolate(fe, V).vector()
    # plot(f, range_min=0., range_max=1.)
    plot(f)

interactive()
