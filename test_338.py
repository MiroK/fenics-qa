from dolfin import *

def foo(mesh):
    x = SpatialCoordinate(mesh)
    f = sin(pi/2*(x[0]+x[1]))
    f_expr = Expression('sin(pi/2*(x[0]+x[1]))', degree=4)

    W = VectorFunctionSpace(mesh, 'CG', 1)
    g0 = grad(f)
    gh = project(g0, W)

    V = FunctionSpace(mesh, 'CG', 1)
    u, v = TrialFunction(V), TestFunction(V)
    a = inner(grad(u), grad(v))*dx
    L = inner(grad(v), gh)*dx

    A, b = assemble_system(a, L)

    null_space_vector = interpolate(Constant(1), V).vector()
    null_space_vector[:] /= null_space_vector.norm('l2')

    null_space = VectorSpaceBasis([null_space_vector])
    as_backend_type(A).set_nullspace(null_space)

    null_space.orthogonalize(b);

    PETScOptions.set("ksp_type", "cg")
    PETScOptions.set("pc_type", "gamg")
    PETScOptions.set("mg_coarse_ksp_type", "preonly")
    PETScOptions.set("mg_coarse_pc_type", "svd")
    PETScOptions.set("ksp_rtol", 1.0e-8)

    # Create Krylov solver and set operator
    solver = PETScKrylovSolver()
    solver.set_operator(A)

    # Set PETSc options on the solver
    solver.set_from_options()

    # Create solution Function
    uh = Function(V)

    # Solve
    niters = solver.solve(uh.vector(), b)

    return mesh.hmin(), errornorm(f_expr, uh, 'H10'), niters

# ----------------------------------------------------------------------------------------

if __name__ == '__main__':
    fmt = 'h=%g, error=%g, rate=%.2f, niters=%d'

    h0, errors = 0., []
    for ncells in (2**i for i in range(3, 10)):
        mesh = UnitSquareMesh(ncells, ncells)

        h, e, niters = foo(mesh)

        if errors:
            print fmt % (h, e, ln(e/errors[-1])/ln(h/h0), niters)
        h0 = h
        errors.append(e)
