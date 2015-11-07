from dolfin import *
from petsc4py import PETSc
import numpy as np

mesh = UnitSquareMesh(4, 4)
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)

u = TrialFunction(V)
v = TestFunction(V)
p = TrialFunction(Q)
q = TestFunction(Q)

a = div(v) * p * dx

# Define boundary condition
u0 = Constant((0.0, 0.0))
bc = DirichletBC(V, u0, 'on_boundary')

M = PETScMatrix()
assemble(a, tensor=M, keep_diagonal=True)

# This operation succeeds
# bc.zero(M)

# Extract the boundary dofs
bcdict = bc.get_boundary_values()
bcinds = np.array(bcdict.keys(), dtype='int32')

# The matrix is rectangular and diagonal that is kep is the smaller square one
N = M.mat()
for row in np.sort(bcinds):
    cols = N.getRow(row)[0]
    print row, row in cols

print bcinds[bcinds < N.size[1]]

N.zeroRowsLocal([0])
