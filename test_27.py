from dolfin import *
from scipy.sparse import csr_matrix
import petsc4py
petsc4py.init()
from petsc4py import PETSc

mesh = UnitIntervalMesh(100)

V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)
bc = DirichletBC(V, Constant(0.), DomainBoundary())

a = inner(grad(u), grad(v))*dx
m = inner(u, v)*dx
L = inner(Constant(0.), v)*dx

A, M = PETScMatrix(), PETScMatrix()
b = PETScVector()
assemble_system(a, L, bc, A_tensor=A, b_tensor=b)
assemble_system(m, L, bc, A_tensor=M, b_tensor=b)

eigensolver = SLEPcEigenSolver(A, M)
eigensolver.parameters['spectrum'] = 'smallest real'
eigensolver.solve(4)

dolfin_eigs = [eigensolver.get_eigenvalue(i)
               for i in range(eigensolver.get_number_converged())]

# Now run the matrices through numpy-csr-pets4py-dolfin chain
A_np = csr_matrix(A.array())
csr = (A_np.indptr, A_np.indices, A_np.data)
A_petsc = PETSc.Mat().createAIJ(size=A_np.shape, csr=csr)
A_petsc.assemble()

M_np = csr_matrix(M.array())
csr = (M_np.indptr, M_np.indices, M_np.data)
M_petsc = PETSc.Mat().createAIJ(size=M_np.shape, csr=csr)
M_petsc.assemble()

AA = PETScMatrix(A_petsc)
MM = PETScMatrix(M_petsc)

eigensolver = SLEPcEigenSolver(AA, MM)
eigensolver.parameters['spectrum'] = 'smallest real'
eigensolver.solve(4)

np_eigs = [eigensolver.get_eigenvalue(i)
           for i in range(eigensolver.get_number_converged())]

for d, n in zip(dolfin_eigs, np_eigs):
    print d, n
