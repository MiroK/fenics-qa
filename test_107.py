from dolfin import *

mesh = UnitSquareMesh(10, 10)
gdim = mesh.geometry().dim()

V = VectorFunctionSpace(mesh, 'CG', 2)
P = FunctionSpace(mesh, 'CG', 1)
M = MixedFunctionSpace([V, P])

u, p = interpolate(Expression(('x[1]', 'x[0]', 'x[0]+x[1]')), M).split()

sigma = -p*Identity(gdim) + 2*sym(grad(u))
n = FacetNormal(mesh)

for i in range(gdim):
    traction_i = (dot(sigma, n))[0]
    ti = project(traction_i, P)
    plot(ti, title=str(i))
interactive()
