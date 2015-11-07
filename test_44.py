from dolfin import *
import numpy as np

mesh = UnitSquareMesh(4, 5)

X = mesh.coordinates().reshape((-1, 2))

sort0 = np.argsort(X[:, 0])
X0 = np.zeros_like(X)
for new, old in enumerate(sort0):
    X0[new, :] = X[old, :]

sort1 = np.argsort(X[:, 1])
Y = np.zeros_like(X0)
for new, old in enumerate(sort1):
    Y[new, :] = X0[old, :]

print Y

# Now Y is a sorted such that Y(i, j) = i*h, j*h
assert all(np.linalg.norm(Y[i] - X[sort0[sort1[i]]]) < 1E-15
           for i in range(len(Y)))

V = FunctionSpace(mesh, 'CG', 1)
u = interpolate(Expression('x[0]'), V)

v2d = vertex_to_dof_map(V)
u_values = u.vector().array()
for i in range(len(X)):
    print Y[i], u_values[v2d[sort0[sort1[i]]]]

plot(u)
interactive()




