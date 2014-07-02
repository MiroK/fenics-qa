from dolfin import *
import numpy as np

mesh = UnitSquareMesh(4, 4)
V = VectorFunctionSpace(mesh, 'CG', 1)

e0 = interpolate(Constant((1, 2)), V)
normalize(e0.vector(), 'l2')

e1 = interpolate(Constant((2, 1)), V)
normalize(e1.vector(), 'l2')

basis_f = [e0, e1]
basis_v = [e.vector() for e in basis_f]

basis = VectorSpaceBasis(basis_v)
print 'Is orthogonal', basis.is_orthogonal()
D = np.zeros((2, 2))
for i, ei in enumerate(basis_v):
    for j, ej in enumerate(basis_v):
        D[i, j] = ei.inner(ej)
print D


print 'Is Orthonormal', basis.is_orthonormal()
e0 = interpolate(Constant((1, 0)), V)

e1 = interpolate(Constant((0, 2)), V)

basis_f = [e0, e1]
basis_v = [e.vector() for e in basis_f]

D = np.zeros((2, 2))
for i, ei in enumerate(basis_v):
    for j, ej in enumerate(basis_v):
        D[i, j] = ei.inner(ej)
print D
