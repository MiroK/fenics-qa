from dolfin import Expression, plot, UnitIntervalMesh
from ufl.algorithms import estimate_total_polynomial_degree as get_degree


class MyExpr(Expression):
    def __init__(self, a, **kwargs):
        self.a = a

    def eval(self, value, x):
        value[0] = x[0] + self.a


f = MyExpr(10, degree=10)
assert get_degree(f) == 10

mesh = UnitIntervalMesh(100)
plot(f, mesh=mesh, interactive=True)
