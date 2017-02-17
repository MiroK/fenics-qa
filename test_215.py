from dolfin import *

parameters['form_compiler']['cpp_optimize'] = True
parameters['form_compiler']['optimize'] = True

mesh = UnitSquareMesh(10, 10)
V = VectorFunctionSpace(mesh, 'CG', 1)
v = interpolate(Expression(('x[0]', '2*x[1]'), degree=1), V)
n = FacetNormal(mesh)

def myjump(v, n, tol=1E-12):
    j = v('+')[0]*n('+')[0] + v('-')[0]*n('-')[0]
    return conditional(abs(n('-')[0]) > tol, j/n('-')[0], j)

print assemble(myjump(v, n)*dS)
