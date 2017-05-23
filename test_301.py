from dolfin import *
import numpy as np

mesh = UnitSquareMesh(10, 10)
V = TensorFunctionSpace(mesh, 'CG', 1)
v = TestFunction(V)

u = Function(V)
U = u.vector()
U.set_local(np.random.rand(U.local_size()))
U.apply('insert')

q, w, e, r = u.split()
u0 = as_matrix(((q, w),
                (e, r)))
assert id(w) == id(u0[0, 1]) 

b = assemble(inner(u0 - u, v)*dx)
print b.norm('linf')

# -----------------------------------

u0 = as_matrix(((q, e),
                (w, r)))

b = assemble(inner(u0 - u, v)*dx)
print b.norm('linf')
