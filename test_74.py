from dolfin import *

mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

a = inner(grad(u), grad(v))*dx + inner(u, v)*dx
m = inner(u, v)*dx

A, M = PETScMatrix(), PETScMatrix()

assemble(a, A)
assemble(m, M)

eigensolver = SLEPcEigenSolver(A,M)
eigensolver.parameters['spectrum'] = 'smallest real'
eigensolver.parameters['tolerance'] = 1.e-14
eigensolver.parameters['problem_type'] = 'gen_hermitian'

eigensolver.solve(5)

for i in range(eigensolver.get_number_converged()):
    r, c, rx, cx = eigensolver.get_eigenpair(i)
    print "%d:  %g" % (i, r) 

