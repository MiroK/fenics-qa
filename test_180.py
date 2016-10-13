import numpy as _np
from scipy.sparse import identity
from petsc4py import PETSc as _PETSc
from dolfin import *

#dummy params
side = int(1e1)
dim = side*side

#rank of present proc and number of cores
comm = _PETSc.COMM_WORLD
rank = comm.Get_rank()
cpu = comm.Get_size()

# matrix definition
A = _PETSc.Mat()
A.create(comm)
A.setSizes([dim,dim])
A.setType('aij')

#dummy matrix
scipyMat = identity(dim)
csr=(scipyMat.indptr,scipyMat.indices,scipyMat.data)
A.setPreallocationCSR(csr)

x,b = A.getVecs()
#define dummy right hand side b
numpRHS = _np.linspace(0.,9.,dim)[b.getOwnershipRange()[0]:b.getOwnershipRange()[1]]
b = _PETSc.Vec().createWithArray(numpRHS,comm=comm)
# set zeros to the solution vector
x = _PETSc.Vec().createWithArray(_np.zeros(dim)[b.getOwnershipRange()[0]:b.getOwnershipRange()[1]],comm=comm) 

SOLVE = True #(or False)
if SOLVE:
 petscSol, petscRHS, petscMat = PETScVector(x), PETScVector(b), PETScMatrix(A)
 comm.barrier() 
 print 'local norm', _np.sum(_np.abs(x.array-b.array))
 solve(petscMat, petscSol, petscRHS, "bicgstab", "icc")
 # verify that x == b
 print 'local norm', _np.sum(_np.abs(x.array-b.array))

print 'Vector', b.getOwnershipRange(), x.getOwnershipRange(), 'Matrix', A.getOwnershipRange(), 'rank', rank
