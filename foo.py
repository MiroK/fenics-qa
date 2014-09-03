from sympy import *
import numpy as np

X = symbols('x y z')
x, y, z = X

def box_integral(f):
    i0 = integrate(f, (x, -2, 2))
    i1 = integrate(i0, (y, -1, 1))
    i2 = integrate(i1, (z, -1, 1))
    return i2

delta = np.eye(3)

inertia_tensor = [[box_integral((x**2 + y**2 + z**2)*delta[i, j] - X[i]*X[j])
                   for j in range(3)] for i in range(3)]

print inertia_tensor
