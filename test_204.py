from dolfin import *

n = 50
mesh = UnitSquareMesh(n, n)
V = FunctionSpace(mesh, "Lagrange", 1)

def boundary(x):
    return x[0] < DOLFIN_EPS or x[1] < DOLFIN_EPS 
u0 = Constant(0.0)
bc = DirichletBC(V, u0, boundary)

r   = Expression("1.", degree=1)

f = Expression("10.", degree=1)

w = TrialFunction(V)
v = TestFunction(V)
a = inner(grad(w), grad(v))*dx
L = f*v*dx + r*v*ds

A = assemble(a)
l = assemble(L)

bc.apply(A)
bc.apply(l)

y = Function(V)
solve(A, y.vector(), l)

file = File('Poisson.pvd')
file << y

