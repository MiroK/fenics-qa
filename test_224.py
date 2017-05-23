from dolfin import *
import numpy as np

mesh = UnitCubeMesh(20, 20, 20)
V = FunctionSpace(mesh, 'CG', 1)

f = interpolate(Expression('2*x[0]+3*x[1]+x[2]', degree=1), V)

Q = VectorFunctionSpace(mesh, 'DG', 0)
q = TestFunction(Q)

MinvL = (1/CellVolume(mesh))*inner(grad(f), q)*dx
x = assemble(MinvL)
grad_f = Function(Q, x)

grad_f0 = interpolate(Constant((2, 3, 1)), Q)

print (grad_f.vector() - grad_f0.vector()).norm('linf')
