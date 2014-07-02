from dolfin import *

mesh = UnitSquareMesh(10, 10)

V = FunctionSpace(mesh, 'CG', 1)

ui = interpolate(Expression('sin(x[0])'), V)
vi = interpolate(Expression('cos(x[1])'), V)
uv_i = interpolate(Expression('sin(x[0]) + cos(x[1])'), V)
uv_p = project(Expression('sin(x[0]) + cos(x[1])'), V)

# Build uvi=u+v by adding expansion coefficients
uvi = Function(V)
uvi.vector().axpy(1, ui.vector())
uvi.vector().axpy(1, vi.vector())
# Check the error
uvi.vector()[:] -= uv_i.vector()
print uvi.vector().norm('linf')

# Build uvp = u+v by projection
uvp = project(ui + vi, V)
# Check the error
uvp.vector()[:] -= uv_p.vector()
print uvp.vector().norm('linf')
