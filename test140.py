from dolfin import *


class SymTransZ(Expression):
    'Given u: (x, y) --> R create v: (x, y, z) --> R, v(x, y, z) = u(x, y).'
    def __init__(self, u):
        self.u = u

    def eval(self, values, x):
        values[0] = self.u(x[0], x[1])

mesh2d = RectangleMesh(-1, -1, 1, 1, 8, 8)
V = FunctionSpace(mesh2d, 'CG', 1)
u = interpolate(Expression('x[0]*x[0] + x[1]*x[1]'), V)
plot(u)

mesh3d = BoxMesh(-1, -1, -1, 1, 1, 1, 16, 16, 5)
W = FunctionSpace(mesh3d, 'CG', 1)
U = SymTransZ(u)
v = interpolate(U, W)

plot(v)
interactive()
