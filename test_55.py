from dolfin import *

mesh = UnitSquareMesh(100, 1000)

V = VectorFunctionSpace(mesh, 'CG', 2)
# Some velocity
v = interpolate(Expression(('sin(pi*(x[0]+x[1]))',
                            'cos(pi*(x[0]-x[1]))')), V)

# Stabilization definition following
# Stability of the SUPG Finite Element Method for Transient Advection-Diffusion Problems
# but ignore polynomial degree
h = CellSize(mesh)
v_mag = sqrt(inner(v, v))
# This is simplified. The criteria should depend on Peclet number. The point is
# that there is no stabilization needed if the flow is mellow. Without the
# conditional we divide almost by zero which makes the solver unhappy
tau = conditional(v_mag > 1E-10, h/2/v_mag, 0)

# Project to DG0, not really needed. The above can be evaluated in any function
# space
S = FunctionSpace(mesh, 'DG', 0)
tau = project(tau, S)

plot(tau)
interactive()
