from ufl.algorithms import estimate_total_polynomial_degree as get_degree
from dolfin import *
mesh = UnitIntervalMesh(2)

V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

k1 = Expression('sin(x[0])*exp(x[0])', degree=2)
a1 = k1*u*v*dx
A1 = as_backend_type(assemble(a1)).array()
print 'a1 degree', get_degree(a1)

k2 = Expression('sin(x[0])*exp(x[0])', degree=3)
a2 = k2*u*v*dx
A2 = as_backend_type(assemble(a2)).array()
print 'a2 degree', get_degree(a2)

b = k2*u*v*dx(metadata={'quadrature_degree': 4})
B = as_backend_type(assemble(b)).array()

print (A2-B)
