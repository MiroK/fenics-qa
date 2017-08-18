from dolfin import *
import numpy as np
import pdb

def LeastSquare(A, b):
    from numpy import linalg
    result = linalg.lstsq(A, b)[0]
    return result

mesh = UnitSquareMesh(5,5)

V = FunctionSpace(mesh, "CG", 1)
W = FunctionSpace(mesh, "CG", 2)

f = Expression("sin(3.14*x[0])", degree=2)
u = TrialFunction(V)
v = TestFunction(V)
w = TestFunction(W)

bc = DirichletBC(W, Constant(1.0), DomainBoundary())
print 'Index and values of DOFs on boundary:', bc.get_boundary_values()

# Forms for linear system, with different test/trial functions
a = u*w*dx + dot(grad(u), grad(w))*dx
L = f*w*dx

# Create linear system (matrices, transpose operator, RHS vector)
A = assemble(a)
b = assemble(L)

Amat = as_backend_type(A).mat()
bvec = as_backend_type(b).vec()

bc.apply(A, b)

# Solve and plot
u = Function(V)
u.vector()[:] = LeastSquare(A.array(), b.array())
plot(u, title='solution')

interactive()
