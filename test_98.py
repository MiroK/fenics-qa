from dolfin import *
import numpy as np

# The requirement on input is that it is a vector valued function in 3d
# For simplicity it is takes as Expression here
class Foo(Expression):
    def eval(self, values, x):
        values[0] = x[0]*x[2]
        values[1] = x[1]*x[2]
        values[2] = 0.
    def value_shape(self):
        return (3, )

f = Foo()
z = 3.

mesh = UnitSquareMesh(10, 10)
V = VectorFunctionSpace(mesh, 'CG', 1)
x = mesh.coordinates().reshape((-1, 2))
v2d = vertex_to_dof_map(V)

values = np.zeros(V.dim(), dtype=float)
for vertex in range(x.shape[0]):
    y = f(x[vertex, 0], x[vertex, 1], z)
    dofs = [v2d[2*vertex+component] for component in range(2)]
    values[dofs] = y

f = Function(V)
f.vector()[:] = values

plot(f)
interactive()

