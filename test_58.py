from dolfin import *

class Foo(Expression):
    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2

    def eval(self, values, x):
        if x[0] < 0.25 + DOLFIN_EPS:
            values[0] = self.p1(x)
        elif x[0] < 0.5 + DOLFIN_EPS:
            values[0] = self.p2(x)
        else:
            values[0] = 0

mesh = UnitSquareMesh(20, 20)

p1 = Expression('x[0]', degree=1, domain=mesh)
p2 = Expression('1-x[0]', degree=1, domain=mesh)
f = Foo(p1, p2)

g = Expression('x[0] < 0.25 + DOLFIN_EPS ? p1 : (x[0] < 0.5 + DOLFIN_EPS ? p2 : 0)',
               p1=p1, p2=p2)

V = FunctionSpace(mesh, 'CG', 1)
f = interpolate(f, V)
g = interpolate(g, V)

plot(f, title='f')
plot(g, title='g')
interactive()
