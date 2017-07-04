from dolfin import *
from petsc4py import PETSc
from scipy.sparse import csr_matrix, coo_matrix
import numpy as np

mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, 'CG', 1)
u, v = TrialFunction(V), TestFunction(V)
M = assemble(inner(u, v)*dx)

Mmat = as_backend_type(M).mat()
# Checkout PETSc.Viewer.createBinary or createASCII with mtx
mat = coo_matrix(csr_matrix(Mmat.getValuesCSR()[::-1]))

np.savetxt('foo.mtx',
           np.c_[mat.row, mat.col, mat.data],
           fmt=['%d', '%d', '%.16f'])
