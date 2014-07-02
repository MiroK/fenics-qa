from dolfin import *

u = VectorConstant(triangle)
u_x, u_y = split()

v = u/sqrt(inner(u, u))

norm = sqrt(u_x**2 + u_y**2)
w = v/norm

mesh = UnitSquareMesh(20, 20)
V = VectorFunctionSpace(mesh, 'CG', 1)
plot(project(u, V), title='u')
plot(project(v, V), title='v')
plot(project(w, V), title='w')
interactive()
