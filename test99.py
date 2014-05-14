from dolfin import *
from numpy import zeros, array

parameters['linear_algebra_backend']= 'uBLAS'

mesh = UnitSquareMesh(10, 10)

V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

f = Expression('sin(pi*x[0])*x[1]')
a = inner(grad(u), grad(v))*dx
L = inner(f, v)*dx
bc = DirichletBC(V, Constant(0.), DomainBoundary())

A, b = assemble_system(a, L, bc)
columns, values = A.getrow(0)
values = zeros(len(values))
A.setrow(0, array(columns, 'uintp'), values)
A.apply('insert')

u = Function(V)
solver = UmfpackLUSolver()
solver.parameters['report'] = True
solver.solve(A, u.vector(), b)
status = solver.umfpack_check_status()
print status
