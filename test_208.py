from dolfin import *
from mshr import *
import sys

# mesh = UnitSquareMesh(*[int(sys.argv[1])]*2)
mesh = generate_mesh(Circle(Point(0, 0), 1), int(sys.argv[1]))

V = VectorFunctionSpace(mesh, 'CG', 1)
v = interpolate(Expression(('x[0]', 'x[1]'), degree=1), V)


Q = FunctionSpace(mesh, 'DG', 0)
q = TestFunction(Q)
n = FacetNormal(mesh)
# Projection step accounting for the mass matrix inverse
# TrialFunction(Q)*q*ds -> FacetArea(mesh)
h = FacetArea(mesh)
v_dot_n = assemble(dot(v, n)*q/h*ds)
# As function
v_dot_n = Function(Q, v_dot_n)

# Check: \int_{\Omega} div(v)*dx = \int_{\partial\Omega} v.n*ds
value0 = assemble(div(v)*dx)
value = assemble(v_dot_n*ds)

info('Error %g' % abs(value0-value))
