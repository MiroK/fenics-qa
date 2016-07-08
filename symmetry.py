from dolfin import *

def eigenvalues(V, bcs, tol, symmetry):
    u = TrialFunction(V)
    v = TestFunction(V)
    dummy = inner(Constant((1, 0)), v)*dx

    # Assemble matrix for LHS
    a = inner(curl(u), curl(v))*dx
    asm = SystemAssembler(a, dummy, bcs)
    A = PETScMatrix(); asm.assemble(A)

    # Assemble matrix for RHS
    b = inner(u, v)*dx
    asm = SystemAssembler(b, dummy)
    B = PETScMatrix(); asm.assemble(B)
    b = PETScVector(); asm.assemble(b)
    
    if symmetry: [bc.zero_columns(B, b) for bc in bcs]
    
    [bc.zero(B) for bc in bcs]

    return bool(True*A.mat().isSymmetric(tol)), bool(True*B.mat().isSymmetric(tol))


def compare_eigenvalues(mesh, tol, symmetry):
    nedelec_V   = FunctionSpace(mesh, "N1curl", 1)
    nedelec_bcs = [DirichletBC(nedelec_V, Constant((0.0, 0.0)), DomainBoundary())]
    nedelec_eig = eigenvalues(nedelec_V, nedelec_bcs, tol, symmetry)

    lagrange_V   = VectorFunctionSpace(mesh, "Lagrange", 1)
    lagrange_bcs = [DirichletBC(lagrange_V.sub(1), 0, "near(x[0], 0) || near(x[0], pi)"),
                    DirichletBC(lagrange_V.sub(0), 0, "near(x[1], 0) || near(x[1], pi)")]
    lagrange_eig = eigenvalues(lagrange_V, lagrange_bcs, tol, symmetry)

    ans = ', '.join(['Ned symmetry A:%s, B:%s' % nedelec_eig,
                     'Lag symmetry A:%s, B:%s' % lagrange_eig])
    return ans

# -----------------------------------------------------------------------------

symmetry = True

mesh = RectangleMesh(Point(0, 0), Point(pi, pi), 40, 40)
print("\ndiagonal mesh")
for tol in [10.**-i for i in range(12)]:
    print '%.2E'% tol, compare_eigenvalues(mesh, tol, symmetry)

mesh = RectangleMesh(Point(0, 0), Point(pi, pi), 40, 40, "crossed")
print("\ncrossed mesh")
for tol in [10.**-i for i in range(12)]:
    print '%.2E' % tol, compare_eigenvalues(mesh, tol, symmetry)
