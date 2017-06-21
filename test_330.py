from dolfin import *

mesh = UnitSquareMesh(10, 10)

V = VectorElement('Lagrange', triangle, 2)
Q = FiniteElement('Lagrange', triangle, 1)
W = MixedElement([V, Q])
W = FunctionSpace(mesh, W)

# V x Q valued expression
w = interpolate(Expression(('x[1]', 'x[0]', 'x[1]+x[0]'), degree=1), W)
u, p = w.split()

# Work with components instead
w0 = Function(W)
u0, p0 = w0.split()

u_init = Expression(('x[1]', 'x[0]'), degree=1)
p_init = Expression('x[1]+x[0]', degree=1)
# to u0 from u_init interpolated in V
assign(u0, interpolate(u_init, W.sub(0).collapse()))  
# to p0 from p_init interpolated in Q
assign(p0, interpolate(p_init, W.sub(1).collapse()))  

# Same result
print assemble(inner(u0-u, u0-u)*dx), assemble(inner(p0-p, p0-p)*dx)
