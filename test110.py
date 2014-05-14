from dolfin import *


class PhiExpression(Expression):
    def __init__(self, F):
        Expression.__init__(self)
        mesh = F.function_space().mesh()
        W = FunctionSpace(mesh, 'DG', 0)
        self.divF = project(div(F), W)
        self.r_mag = \
            Expression('sqrt((x[0]-a)*(x[0]-a)+(x[1]-b)*(x[1]-b))', a=0, b=0)

    def eval(self, values, x):
        self.r_mag.a = x[0]
        self.r_mag.b = x[1]
        values[0] = assemble(self.divF/self.r_mag*dx)

mesh = UnitSquareMesh(20, 20)
V = FunctionSpace(mesh, 'CG', 1)
M = VectorFunctionSpace(mesh, 'CG', 1)
F = interpolate(Expression(('x[0]', 'x[0]*x[1]')), M)

phi = interpolate(PhiExpression(F), V)
plot(F, mesh=mesh)
plot(phi)
interactive()
