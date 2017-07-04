from dolfin import *

mesh = UnitIntervalMesh(30)
# A real space in fenics represents a constant function. A space of functions 
# {f0, f1, f2} can be represented in fenics using f0*b0 \times f1*b1 \times f2*b2
# where bi are the basis functions of vector space of constants 

# Function f0 must be something that fenics understands, here I go for expressions
f0 = Expression('sin(k*x[0])', degree=6, k=pi)      # degree for getting quadrature 
f1 = Expression('sin(k*x[0])', degree=6, k=2*pi)
f2 = Expression('sin(k*x[0])', degree=6, k=3*pi)

fs = [f0, f1, f2]

# The vector space of constants
V = VectorFunctionSpace(mesh, 'R', 0, dim=len(fs))
us = TrialFunction(V)
vs = TestFunction(V)
# f*u
u = sum(fi*bi for fi, bi in zip(fs, us))  # custom trial functions
# f*v
v = sum(fi*bi for fi, bi in zip(fs, vs))  # custom test functions

# You can use it in a form now

a = inner(u, v)*dx
print assemble(a).array()



