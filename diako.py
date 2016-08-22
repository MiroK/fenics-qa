from dolfin import *
import numpy as np

# Source term
f0 = Expression('sin(k*x[0])', k=2*pi, degree=4)
f1 = Expression('sin(k*x[0])', k=1*pi, degree=4)

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

# Create mesh and finite element
for ncells in [32, 64, 128, 256, 512]:
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
    solver = PETScKrylovSolver('cg')
    solver.set_operator(A)
    # Create vector that spans the null space and normalize
    null_vec = interpolate(Constant(1), V).vector()
    null_vec *= 1.0/null_vec.norm('l2')
    # Create null space basis object and attach to PETSc matrix
    null_space = VectorSpaceBasis([null_vec])
    as_backend_type(A).set_nullspace(null_space)
    # Orthogonalize RHS vector b with respect to the nullspace
    alpha = null_vec.inner(b)
    b.axpy(-alpha, null_vec)
    # Solve
    solver.solve(uh.vector(), b)

    # # Stage two: neumann
    # V = FunctionSpace(mesh, 'CG', 1)
    # u = TrialFunction(V)
    # v = TestFunction(V)

    # a = inner(grad(u), grad(v))*dx
    # A = assemble_system(a)
    # # Rhs
    # b = interpolate(Constant(1), V).vector()
    # b *= 1.0/b.norm('l2')
    # # null_vec

    print 'h = %4f e = %6f' % (mesh.hmin(), errornorm(u_exact, uh, 'H1'))

    #x0[:] = b1

# Visual
plot(uh, title='numeric')
plot(u_exact, mesh=mesh, title='exact')

e = interpolate(u_exact, V)
e.vector().axpy(-1, uh.vector())
e.vector()[:] += e.vector().min()

Ae = e.vector().copy()
A.mult(e.vector(), Ae)
print Ae.norm('linf'), null_vec.inner(Ae)

plot(e, title='error')

interactive()
