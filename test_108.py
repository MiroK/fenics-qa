from dolfin import *
import numpy as np

A = np.array([[10, 2],
              [1, -3]])
Sym = np.array([[1, 0.5],
                [0.5, 2]])

mesh = UnitTriangleMesh()
V = TensorFunctionSpace(mesh, 'CG', 1)
S = TensorFunctionSpace(mesh, 'CG', 1, symmetry=True)

# See if T can be represented exactly in S
test = lambda T, S: np.allclose(interpolate(Constant(T), S)(0.5, 0.5).reshape((2, 2))-T,
                                1E-13)
# Sanity
assert test(A, V)
assert test(Sym, V)
# Correcly comes out as False - can't represent not symm tensor. 
print test(A, S)           
print test(Sym, S)         # BUG: False despite symmetry. True only if Sym[1, 1]
Sym[1, 1] = 0
print test(Sym, S)         
