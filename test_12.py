from dolfin import *
import time


for N in (128, 192, 256, 384):
    mesh = UnitSquareMesh(N, N)

    V = FunctionSpace(mesh, 'CG', 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    bc = DirichletBC(V, Constant(0), DomainBoundary())

    a = inner(grad(u), grad(v))*dx
    L = inner(Constant(1), v)*dx

    u = Function(V)
    problem = LinearVariationalProblem(a, L, u, bc)
    solver = LinearVariationalSolver(problem)
    solver.linear_solver = 'petsc'
    solver.parameters['symmetric'] = True
    
    start = time.time()
    solver.solve()
    print time.time() - start

