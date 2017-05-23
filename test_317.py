from dolfin import *
import matplotlib.pyplot as plt
import numpy as np

mesh = UnitIntervalMesh(10)

V = FunctionSpace(mesh, 'CG', 2)
f = interpolate(Expression('x[0]*x[0]', degree=2), V)

y = f.vector().array()

plt.figure()
plt.plot(y, linestyle='none', marker='o')  

# As you noticed the coefficients are not ordered according 
# to x coord of dofs. Get the x axis 'right' by e.g.
x = interpolate(Expression('x[0]', degree=1), V).vector().array()
# or
x0 = V.tabulate_dof_coordinates()
assert np.linalg.norm(x-x0) < DOLFIN_EPS

plt.figure()
plt.plot(x, y, linestyle='none', marker='o')

# This is better, but what if you wanted a line plot
plt.figure()
plt.plot(x, y)
# You see that the plot isn't quite right; the interior dofs of CG2 element are an 
# issue. So finally
indices = np.argsort(x)
x = x[indices]
y = y[indices]
plt.figure()
plt.plot(x, y)

plt.show()
