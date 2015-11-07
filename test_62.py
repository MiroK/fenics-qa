from dolfin import *

mesh = UnitSquareMesh(10, 10)
V = FunctionSpace ( mesh, "CG", 1 )

class MyExpression(Expression):
    def eval(self, value, x):
        if x[0] <= .5:
            value[0] = 1.0
        else:
            value[0] = 0.0

    def value_shape(self):
        return ()

u_init = MyExpression()
u = interpolate ( u_init, V )
print u.vector().array()

print interpolate(Expression("x[0] < 0.5 + DOLFIN_EPS ? 1: 0"),
        V).vector().array()
