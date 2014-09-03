from dolfin import *

mesh = UnitSquareMesh(10, 10)

V = VectorFunctionSpace(mesh, 'CG', 1)
Q = FunctionSpace(mesh, 'CG', 1)
M = MixedFunctionSpace([V, Q])

m_exact = Expression(('x[0]', 'x[1]', 'x[0]+x[1]'))
u_exact = Expression(('x[0]', 'x[1]'))
p_exact = Expression('x[0]+x[1]')

# Hypotetical solution of the mixed problem
mh = interpolate(m_exact, M)
# Perturb the solution
mh.vector()[:] += 1
uh, ph = mh.split()

m_error = sqrt(assemble(inner(m_exact - mh, m_exact - mh)*dx))
u_error = errornorm(u_exact, uh, 'l2')
p_error = errornorm(p_exact, ph, 'l2')

print m_error, '%.2e' % abs(m_error - sqrt(3))
print u_error, '%.2e' % abs(u_error - sqrt(2))
print p_error, '%.2e' % abs(p_error - 1)
