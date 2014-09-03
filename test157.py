from dolfin import *
import numpy as np

mesh = UnitSquareMesh(2, 2)

V = VectorFunctionSpace(mesh, 'CG', 1)
u = interpolate(Expression(('x[0]', 'x[1]')), V)

T = TensorFunctionSpace(mesh, 'CG', 1)
t = project(grad(u), T)

print t.vector().norm('l2')

t_array = t.vector().array()
t_array += np.random.rand(len(t_array))
t.vector()[:] = t_array

print t.vector().norm('l2')
