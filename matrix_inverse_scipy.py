from dolfin import *
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import inv
import matplotlib.pyplot as plt

# Choose backend that allows .data()
parameters.linear_algebra_backend = 'uBLAS'

mesh = UnitSquareMesh(30, 30)
V = FunctionSpace(mesh, 'CG', 1)

u = TrialFunction(V)
v = TestFunction(V)
a = inner(grad(u), grad(v))*dx
m = inner(u, v)*dx
L = Constant(0.)*v*dx

bc = DirichletBC(V, Constant(0), DomainBoundary())
A, _ = assemble_system(a, L, bc)
M, _ = assemble_system(m, L, bc)

# Prepare A and M for scipy
rows, cols, values = A.data()
A = csr_matrix((values, cols, rows))

rows, cols, values = M.data()
M = csr_matrix((values, cols, rows))

print 'Inverting %d by %d matrix' % (V.dim(), V.dim())
timer = Timer('inverse')
timer.start()
T = inv(M).dot(A)
timer.stop()
print '\t\t done in %g' % timing('inverse')

# Plot
plt.figure()
plt.subplot(131)
plt.spy(A)
plt.subplot(132)
plt.spy(M)
plt.subplot(133)
plt.spy(T)
plt.show()
