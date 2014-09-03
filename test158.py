from dolfin import *

mesh = UnitIntervalMesh(20)

V = FunctionSpace(mesh, 'CG', 1)
M = MixedFunctionSpace([V, V])

# Some function in M to illustrate that components
# will change by assign
m = interpolate(Expression(('x[0]', '1-x[0]')), M)
m0, m1 = split(m)

# Plot for comparison
plot(m0, title='m0 before')
plot(m1, title='m1 before')

# Functions for components
v0 = interpolate(Expression('cos(pi*x[0])'), V)
v1 = interpolate(Expression('sin(pi*x[0])'), V)

# Assign the components
assign(m.sub(0), v0)
assign(m.sub(1), v1)

# See if it worked
plot(m0, title='m0 after')
plot(m1, title='m1 after')
interactive()
