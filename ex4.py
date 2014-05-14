import sys, slepc4py
slepc4py.init(sys.argv)

from slepc4py import SLEPc

from dolfin import *

mesh = UnitSquareMesh(20, 20)

CG = FunctionSpace(mesh, 'CG', 1)
CR = FunctionSpace(mesh, 'CR', 1)
cg = TrialFunction(CR)
cr = TestFunction(CG)

# A is a matrix from CR to CG, m=CG.dim() x n=CR.dim()
a = inner(cg, cr)*dx
A = PETScMatrix()
A = assemble(a, tensor=A)
# print 'A is %dx%d' % (A.size(0), A.size(1))
# print 'A is %dx%d' % (CG.dim(), CR.dim())
A_mat = A.mat()

S = SLEPc.SVD()
S.create()
S.setOperator(A_mat)
S.setType(S.Type.LANCZOS)
S.solve()

nconv = S.getConverged()
if nconv > 0:
    v, u = A_mat.getVecs()    # A*v = u
    # print 'm=%d' % u.size, 'n=%d' % v.size
    print 'Singular values:'
    for i in range(nconv):
        sigma = S.getSingularTriplet(i, u, v)
        print sigma
