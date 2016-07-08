'''
This demo program solves the curl-curl eigenproblem: 

find u in H_0(curl) and \lambda \in \mathbb{R} such that

\int_\Omega curl(u) curl(v) = \lambda \int_\Omega u \cdot v

for all v in H_0(curl), computing the first positive
eigenvalues.

The problem is discretised with Nedelec and piecewise linear
finite elements on a uniform diagonal mesh, and a uniform
crossed mesh. The Nedelec solution approximates the
correct eigenvalues correctly in both cases. The Lagrange
solution gives completely wrong eigenvalues on the diagonal
mesh, and has a spurious eigenvalue with value near 6 for the
crossed mesh.

For more details, see

http://umn.edu/~arnold/papers/icm2002.pdf

Contributed by Patrick E. Farrell and Douglas N. Arnold.
'''

from dolfin import *
import numpy as np

def eigenvalues(V, bcs):

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
    [bc.zero_columns(B, b) for bc in bcs]
    [bc.zero(B) for bc in bcs]

    print 'A symmetry', A.mat().isSymmetric(1E-8)
    print 'B symmetry', B.mat().isSymmetric(1E-8)


    solver = SLEPcEigenSolver(A, B)
    solver.parameters["solver"] = "krylov-schur"
    solver.parameters["spectrum"] = "target magnitude"
    solver.parameters["problem_type"] = "gen_hermitian"
    solver.parameters["spectral_transform"] = "shift-and-invert"
    solver.parameters["spectral_shift"] = 5.5
    neigs = 12
    solver.solve(neigs)

    computed_eigenvalues = []

    for i in range(min(neigs, solver.get_number_converged())):
        (r, _, rv, _) = solver.get_eigenpair(i)
        computed_eigenvalues.append(r)
    return np.sort(np.array(computed_eigenvalues))

def compare_eigenvalues(mesh):
    
    nedelec_V   = FunctionSpace(mesh, "N1curl", 1)
    nedelec_bcs = [DirichletBC(nedelec_V, Constant((0.0, 0.0)), DomainBoundary())]
    nedelec_eig = eigenvalues(nedelec_V, nedelec_bcs)

    lagrange_V   = VectorFunctionSpace(mesh, "Lagrange", 1)
    lagrange_bcs = [DirichletBC(lagrange_V.sub(1), 0, "near(x[0], 0) || near(x[0], pi)"),
                    DirichletBC(lagrange_V.sub(0), 0, "near(x[1], 0) || near(x[1], pi)")]
    lagrange_eig = eigenvalues(lagrange_V, lagrange_bcs)

    true_eig = np.sort(np.array([float(m**2 + n**2) for m in range(6) for n in range(6)]))[1:13]

    np.set_printoptions(formatter={'float': '{:5.2f}'.format})
    print "Nedelec:  ",
    print nedelec_eig 
    print "Lagrange: ",
    print lagrange_eig
    print "Exact:    ",
    print true_eig

mesh = RectangleMesh(Point(0, 0), Point(pi, pi), 40, 40)
print("\ndiagonal mesh")
x = compare_eigenvalues(mesh)

mesh = RectangleMesh(Point(0, 0), Point(pi, pi), 40, 40, "crossed")
print("\ncrossed mesh")
compare_eigenvalues(mesh)
