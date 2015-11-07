from dolfin import *
import numpy as np

mesh = UnitSquareMesh(3, 3)

V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

bc = DirichletBC(V, Constant(0), DomainBoundary())

a = inner(grad(u), grad(v))*dx
m = inner(u, v)*dx
L = inner(Constant(1), v)*dx

A, M = PETScMatrix(), PETScMatrix()
b = PETScVector()

# See that assemble + bc.apply does not imply symmetry
assemble(a, A)
bc.apply(A)
print 'Is symmetric?', np.linalg.norm(A.array() - A.array().T) < 1E-10

# But assemble system does
assemble_system(a, L, bc, A_tensor=A, b_tensor=b)
assemble_system(m, L, bc, A_tensor=M, b_tensor=b)

print 'Is symmetric?', np.linalg.norm(A.array() - A.array().T) < 1E-10
print 'Is symmetric?', np.linalg.norm(M.array() - M.array().T) < 1E-10

from slepc4py import SLEPc

E = SLEPc.PEP().create()
E.setOperators([A.mat(), M.mat(), A.mat()])
