from dolfin import *
import numpy as np

mesh = UnitSquareMesh(20,20)

parameters["linear_algebra_backend"] = "PETSc"

P2 = VectorFunctionSpace(mesh, "CG", 2)
P1 = FunctionSpace(mesh, "CG", 1)
W = MixedFunctionSpace([P2,P1])

bc = DirichletBC(W.sub(0), Constant((0,0)), DomainBoundary())

(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)

a = inner(grad(u), grad(v))*dx 
L = inner(Constant((0,0)), v)*dx

A = PETScMatrix()

assemble(a, A)
print 'Is symmetric', np.linalg.norm(A.array() - A.array().T) < 1E-11
temp=A.mat()
print 'Is symmetric', temp.isSymmetric(1E-11)
