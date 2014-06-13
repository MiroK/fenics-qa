from dolfin import *
import numpy as np

def is_inside_sphere(x, x0, r):
    'Check if x is inside sphere of radius r centered at x0.'
    return np.sum((x-x0)**2) < r + DOLFIN_EPS

mesh = UnitSquareMesh(20, 20)
x0 = np.array([0.5, 0.5])
r = 0.25

class MyExpression(Expression):
    def eval(self, values, x):
        values[0] = 1 if is_inside_sphere(x, x0, r) else 2

plot(MyExpression(), mesh=mesh)
interactive()


