from dolfin import *
mesh = UnitSquareMesh(25, 25)
V = FunctionSpace(mesh, "CG", 1)
u = TrialFunction(V)
v = TestFunction(V)
a = dot(grad(u), grad(v))*dx
m = u*v*dx

print V.dim()

# Assemble the stiffness matrix and the mass matrix.
A = assemble(a)
M = assemble(m)
A = as_backend_type(A).mat()
M = as_backend_type(M).mat()

from slepc4py import SLEPc
from petsc4py import PETSc

# Setup the eigensolver
solver = SLEPc.EPS().create()
solver.setOperators(A, M)
solver.setType(solver.Type.KRYLOVSCHUR)
solver.setDimensions(3, PETSc.DECIDE)
solver.setWhichEigenpairs(solver.Which.SMALLEST_REAL)
solver.setTolerances(1E-8, 1000)

solver.solve()

nconv = solver.getConverged()
niters = solver.getIterationNumber()
if MPI.rank(mesh.mpi_comm()) == 0:
    print "Number of eigenvalues successfully computed: ", nconv
    print "Iterations used", niters

eigenvalues = []
for i in range(nconv):
    r = solver.getEigenvalue(i).real  # real part of eigenvalue
    eigenvalues.append(r)

if MPI.rank(mesh.mpi_comm()) == 0:
    print "Smallest positive eigenvalues computed and exact: "
    for i in range(min(nconv, 3)):
        print eigenvalues[i]
