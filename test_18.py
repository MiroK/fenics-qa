from dolfin import *

mesh = UnitSquareMesh(10, 10)

V = FunctionSpace(mesh, 'CG', 1)
v = TestFunction(V)

n = FacetNormal(mesh)
q_vector = Expression(('x[0]', 'x[1]'))
q = conditional(gt(n[1], 0), dot(q_vector, n), Constant(0.0))
# q = conditional(gt(dot(n, Constant((1, 0))), 0), dot(n, q_vector), Constant(0))
L = q*v*ds
