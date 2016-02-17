from dolfin import *

mesh = UnitSquareMesh(10, 10)

V = VectorFunctionSpace(mesh, 'CG', 1)
Q = FunctionSpace(mesh, 'CG', 1)
W = MixedFunctionSpace([V, Q])
B = FunctionSpace(mesh, 'BDM', 1)

fs = map(Function, (V, Q, W, B))

for f in fs:
    print f, f.function_space().ufl_element().family() == 'Mixed'
