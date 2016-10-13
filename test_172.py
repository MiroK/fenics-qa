from dolfin import *

mesh2d = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh2d, 'CG', 1)
f = interpolate(Expression('sin(pi*(x[0]+x[1]))', degree=4), V)

class Foo(Expression):
    def __init__(self, f): self.f = f

    def value_shape(self): return (3, )

    def eval(self, values, x):
        values[0] = 0.0
        values[1] = 0.0
        values[2] = 0.0 if x[2] < 0.5 else self.f(x[:2])

mesh3d = UnitCubeMesh(10, 10, 10)
V = VectorFunctionSpace(mesh3d, 'CG', 1)
v = interpolate(Foo(f), V)
plot(v)
interactive()

