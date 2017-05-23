from dolfin import *
from petsc4py import PETSc

def mat_inner(A, B):
    A = A.mat()
    B = B.mat()
    
    comm = B.comm
    BtA = PETSc.Mat(comm)
    B.matTransposeMult(A, BtA)

    d = PETSc.Vec(comm)
    BtA.getDiagonal(d)
    return d.sum()

if __name__ == '__main__':
    mesh = UnitSquareMesh(100, 100)
    V = FunctionSpace(mesh, 'CG', 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    A, M = PETScMatrix(), PETScMatrix()
    assemble(inner(grad(u), grad(v))*dx, A)
    assemble(inner(u, v)*dx, M)

    i = mat_inner(A, M)
    print i
