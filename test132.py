from dolfin import *

mesh = UnitSquareMesh(100, 100)

V = FunctionSpace(mesh, 'CG', 1)

u = interpolate(Constant(3), V)
v = interpolate(Constant(4), V)
w = Function(V)

U = u.vector()
V = v.vector()
W = w.vector()

W.axpy(1, U)
W *= V

plot(w)
interactive()
