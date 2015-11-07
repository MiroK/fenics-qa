#
# - laplace(u) = f in [0, 1]^2
#           u  = 0 on the boundary
#         u(p) = 0 for some points in the domain
#
#  Manually assemble the saddle point system
#  [[A, B], [B.T, 0]]  where B are the points constains
#

from dolfin import *
from pyamg.krylov import gmres
from pyamg import smoothed_aggregation_solver
from scipy.sparse import csr_matrix, bmat, eye
from scipy.sparse.linalg import spsolve, LinearOperator
from scipy.linalg import eigvalsh
import numpy as np

        
def mat_to_csr(mat):
    '''Convert any dolfin.Matrix to csr matrix in scipy.'''
    assert MPI.size(mpi_comm_world()) == 1, 'mat_to_csr assumes single process'
    rows = [0]
    cols = []
    values = []
    for row in range(mat.size(0)):
        cols_, values_ = mat.getrow(row)
        rows.append(len(cols_)+rows[-1])
        cols.extend(cols_)
        values.extend(values_)

    shape = mat.size(0), mat.size(1)
    
    return csr_matrix((np.asarray(values),
                       np.asarray(cols, dtype='int'),
                       np.asarray(rows, dtype='int')),
                       shape)
import sys
N = int(sys.argv[1])
mesh = UnitSquareMesh(N, N)

V = FunctionSpace(mesh, 'CG', 1)
bc = DirichletBC(V, Constant(0), DomainBoundary())

# A block and the force
u, v = TrialFunction(V), TestFunction(V)
a = inner(grad(u), grad(v))*dx
L = inner(Constant(1), v)*dx
A, b = assemble_system(a, L, bc)
A_ = mat_to_csr(A)

# B block, one column per point. Represents the action of dela @ point
column = assemble(inner(Constant(0), v)*dx)

points = [Point(0.44, 0.21), Point(0.12, 0.4)]
B_ = np.zeros((A_.shape[0], len(points)))
for col, point in enumerate(points):
    delta = PointSource(V, point)
    delta.apply(column)
    B_[:, col] = column.array()
    column.zero()

# The saddle system
AA = bmat([[A_, B_], [B_.T, None]])
bb = np.r_[b.array(), np.zeros(len(points))]

# Preconditioner diag(ML(A), I)
ml = smoothed_aggregation_solver(A_)
PA = ml.aspreconditioner()

def Pmat_vec(x):
    yA = PA.matvec(x[:V.dim()])
    # Identity on part of x from point sources
    return np.r_[yA, x[V.dim():]]

P = LinearOperator(shape=AA.shape, matvec=Pmat_vec)
    
residuals = []
# Run the iterations
x, failed = gmres(AA, bb, M=P, 
                  tol=1e-10, maxiter=1000, residuals=residuals)
uh = Function(V)
uh.vector()[:] = x[:V.dim()]

# Solve
# xx = spsolve(AA, bb)

# Get the coefs of !multipliers
# U = xx[:V.dim()]
# uh = Function(V)
# uh.vector()[:] = U

print len(residuals)-1, residuals[-1]

plot(uh)
interactive()
