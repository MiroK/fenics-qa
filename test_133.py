from dolfin import *

mesh = UnitSquareMesh(10, 10)
V = VectorFunctionSpace(mesh, 'CG', 2)
v = interpolate(Expression(('x[0]', '-x[1]'), degree=1), V)

print [assemble(v[i]*dx) for i in range(2)]
print [[assemble(Dx(v[i], j)*dx) for j in range(2)] for i in range(2)]
