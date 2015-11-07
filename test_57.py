from dolfin import *
from scipy.sparse import csr_matrix
from scipy.io import savemat

# Convert PETScMatrix to csr_matrix
to_csr = lambda matrix: csr_matrix(tuple(as_backend_type(matrix).mat().getValuesCSR()[::-1]),
                                   shape=(A.size(0), A.size(1)))

parameters['linear_algebra_backend'] = 'PETSc'

mesh = UnitSquareMesh(3, 3)
V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)
bc = DirichletBC(V, Constant(0.), 'on_boundary')

a = inner(grad(u), grad(v))*dx
m = inner(u, v)*dx
L = inner(Constant(1), v)*dx

A, b = assemble_system(a, L, bc)
M, _ = assemble_system(m, L, bc)

savemat('A.mat', dict(A=to_csr(A)))
savemat('M.mat', dict(M=to_csr(M)))
savemat('b.mat', dict(b=b.array()))
