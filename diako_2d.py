from dolfin import *
import numpy as np

u_exact = Expression('sin(k0*x[0])+sin(k1*x[1])+C*cos(2*pi*x[0])',
                     k0=1*pi, k1=3*pi, degree=4, C=1.)

f = Expression('k0*k0*sin(k0*x[0])+k1*k1*sin(k1*x[1])+C*4*pi*pi*cos(2*pi*x[0])',
               k0=u_exact.k0, k1=u_exact.k1, C=u_exact.C, degree=4)

# Sub domain for Periodic boundary condition
class PeriodicBoundary(SubDomain):

    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        # return True if on left or bottom boundary AND NOT on one of the two corners (0, 1) and (1, 0)
        return bool((near(x[0], 0) or near(x[1], 0)) and 
                (not ((near(x[0], 0) and near(x[1], 1)) or 
                        (near(x[0], 1) and near(x[1], 0)))) and on_boundary)

    def map(self, x, y):
        if near(x[0], 1) and near(x[1], 1):
            y[0] = x[0] - 1.
            y[1] = x[1] - 1.
        elif near(x[0], 1):
            y[0] = x[0] - 1.
            y[1] = x[1]
        else:   # near(x[1], 1)
            y[0] = x[0]
            y[1] = x[1] - 1.

# Marking edges
left = CompiledSubDomain('near(x[0], 0.)')
right = CompiledSubDomain('near(x[0], 1.)')
bottom = CompiledSubDomain('near(x[1], 0.)')
top = CompiledSubDomain('near(x[1], 1.)')
boundary_domains = [left, right, bottom, top]

# ----------------------------------------------------------------------------

e0, h0 = -1, -1
for ncells in [8, 16, 32, 64, 128]:
    mesh = UnitSquareMesh(ncells, ncells)

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

    boundaries = EdgeFunction('size_t', mesh, 0)
    [domain.mark(boundaries, tag) for tag, domain in enumerate(boundary_domains, 1)]
    for tag in range(1, 5):
        print '\t', assemble(uh*ds(domain=mesh, subdomain_id=tag, subdomain_data=boundaries))
    print '\t\t', assemble(uh*ds)
    
    W = FunctionSpace(mesh, 'CG', 1)
    # Stage two: Dirichlet correction if the orig right-hand side requires it
    # FIXME: THINK ABOUT THIS.
    # FIXME: TensorProduct...
    if abs(alpha) > 1E-13:
        V1 = FunctionSpace(UnitIntervalMesh(ncells), 'CG', 1)
        u = TrialFunction(V1)
        v = TestFunction(V1)

        a = inner(grad(u), grad(v))*dx
        L = inner(Constant(alpha*nnorm), v)*dx

        points = ([0., 0.5], [0.5, 0.])
        comps = []
        for p in points:
            beta = -uh(*p)
            bc = DirichletBC(V1, Constant(beta), 'on_boundary')
            B, b = assemble_system(a, L, bc)
            vh = Function(V1)

            solver = PETScKrylovSolver('cg', 'hypre_amg')
            solver.set_operator(B)
            solver.solve(vh.vector(), b)
            comps.append(vh)
        # Try to make the tensor product
        dofs = W.tabulate_dof_coordinates().reshape((-1, 2))
        U0 = np.array([comps[0](x) for x in dofs[:, 0]])
        U1 = np.array([comps[1](y) for y in dofs[:, 1]])
        wh = Function(W)
        wh.vector()[:] = U0*U1

        uh = interpolate(uh, W)
        # uh.vector().axpy(1., wh.vector())
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

# Error
e = interpolate(u_exact, W)
e.vector().axpy(-1, uh.vector())
plot(e, title='error')

plot(project(grad(e)[0]**2, W), title='grad_x')
plot(project(grad(e)[1]**2, W), title='grad_y')

plot(wh, title='wh')

interactive()
