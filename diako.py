from dolfin import *
import numpy as np

# Source term
f0 = Expression('sin(k*x[0])', k=3*pi, degree=4)
f1 = Expression('sin(k*x[0])', k=7*pi, degree=4)

u_exact = Expression('f0/k0/k0 + f1/k1/k1', f0=f0, k0=f0.k, f1=f1, k1=f1.k,
                     degree=4)
f = f0 + f1# + Constant(1)

class PeriodicBoundary(SubDomain):
    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        return bool(x[0] < DOLFIN_EPS and x[0] > -DOLFIN_EPS and on_boundary)

    # Map right boundary (H) to left boundary (G)
    def map(self, x, y):
        if near(x[0], 1):
            y[0] = x[0] - 1.0

# ----------------------------------------------------------------------------

e0, h0 = -1, -1
for ncells in [32, 64, 128, 256, 512, 1024]:
    mesh = UnitIntervalMesh(ncells)

    # Stage one: periodic
    V = FunctionSpace(mesh, 'CG', 1, constrained_domain=PeriodicBoundary())
    u = TrialFunction(V)
    v = TestFunction(V)

    a = inner(grad(u), grad(v))*dx
    L = inner(f, v)*dx
    A, b = assemble_system(a, L)
    uh = Function(V)

    # In the first stage we want to avoid constant nullspace
    solver = PETScKrylovSolver('cg', 'hypre_amg')
    solver.set_operator(A)
    # Create vector that spans the null space and normalize
    null_vec = interpolate(Constant(1), V).vector()
    nnorm = null_vec.norm('l2')
    null_vec *= 1.0/nnorm
    # Create null space basis object and attach to PETSc matrix
    null_space = VectorSpaceBasis([null_vec])
    as_backend_type(A).set_nullspace(null_space)
    # Orthogonalize RHS vector b with respect to the nullspace
    alpha = null_vec.inner(b)
    b.axpy(-alpha, null_vec)
    # Solve
    solver.solve(uh.vector(), b)

    beta = uh(0.)
    W = FunctionSpace(mesh, 'CG', 1)
    # Stage two: Dirichlet correction if the orig right-hand side requires it
    if abs(alpha) > 1E-13:
        w = TrialFunction(W)
        v = TestFunction(W)
        bc = DirichletBC(W, Constant(-beta), 'on_boundary')

        a = inner(grad(w), grad(v))*dx
        L = inner(Constant(alpha*nnorm), v)*dx
        B, b = assemble_system(a, L, bc)
        wh = Function(W)

        solver = PETScKrylovSolver('cg', 'hypre_amg')
        solver.set_operator(B)
        solver.solve(wh.vector(), b)

        uh = interpolate(uh, W)
        uh.vector().axpy(1., wh.vector())
    else:
        uh = interpolate(uh, W)
    
    h = mesh.hmin()
    e = errornorm(u_exact, uh, 'H1')
    if e0 > 0:
        rate = ln(e/e0)/ln(h/h0)
        print 'h = %4f e = %6f r = %.2f' % (h, e, rate)
    h0, e0 = mesh.hmin(), e

# Visual
plot(uh, title='numeric')
plot(u_exact, mesh=mesh, title='exact')
interactive()
