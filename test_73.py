from dolfin import *

mesh = UnitSquareMesh(2, 2)

V1= VectorFunctionSpace(mesh, 'CG', 1)
V2 = VectorFunctionSpace(mesh, 'CG', 1)
Q1 = FunctionSpace(mesh, 'CG', 1)
Q2 = FunctionSpace(mesh, 'CG', 1)
VQ = MixedFunctionSpace([V1,V2, Q1, Q2])

class InitialConditions(Expression):
    def eval(self, values, x):
        values[0] = x[1]
        values[1] = x[1]
        values[2] = x[1]
        values[3] = x[1]
        values[4] = 0.01
        values[5] = 0

    def value_shape(self):
        return (6,)

u = InitialConditions()
ic = interpolate(u, VQ)
