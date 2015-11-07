from dolfin import *
import numpy as np
import sys

N = 2**6
mesh = UnitCubeMesh(N, N, N)
V = VectorFunctionSpace(mesh, 'CG', 2)
u = TrialFunction(V)
v = TestFunction(V)

A = PETScMatrix()
assemble(inner(grad(u), grad(v))*dx, tensor=A)

row, col, val = A.mat().getValuesCSR()

print np.iinfo(np.int32).max, np.iinfo(np.int64).max
print len(row), len(col), len(val)
print A.size(0)
