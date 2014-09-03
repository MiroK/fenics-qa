from dolfin import *
import numpy as np

def foo(N):
    mesh = UnitSquareMesh(N, N)

    V = FunctionSpace(mesh, 'CG', 1)

    u = TrialFunction(V)
    v = TestFunction(V)
    bc = DirichletBC(V, Constant(0), DomainBoundary())

    A = PETScMatrix()
    b = PETScVector()

    assemble(inner(grad(u), grad(v))*dx, tensor=A)
    assemble(inner(Constant(0), v)*dx, tensor=b)
    bc.apply(A)

    u = Function(V)
    x = u.vector()

    solver = PETScKrylovSolver('gmres')
    solver.set_operator(A)
    solver.parameters['gmres']['restart'] = 3000
    solver.parameters['monitor_convergence'] = True

    for i in range(1):
        print i,
        b[:] = np.random.rand(V.dim()) - np.random.rand(V.dim())
        bc.apply(b)
        solver.solve(x, b)
        X = np.sort(x.array())
        print X[-1], X[0]

print foo(40)
