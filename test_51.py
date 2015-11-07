from dolfin import *
import petsc4py
petsc4py.init()
from petsc4py import PETSc
import numpy as np

mesh = UnitSquareMesh(20, 20)

V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)
bc = DirichletBC(V, Constant(0.), DomainBoundary())
# lhs
A = PETScMatrix()
assemble(inner(grad(u), grad(v))*dx, tensor=A)
bc.apply(A)
# rhs
B = PETScMatrix()
assemble(inner(u, v)*dx, tensor=B)
# Solver
solver = LUSolver(A)
solver.parameters['reuse_factorization'] = True
# Vectors for solution x, i-th rhs b
n = V.dim()
x, b, ei = [Vector(None, n) for i in range(3)]
# Matrix of solutions x as columns
X = PETSc.Mat().createDense([n, n])
X.setUp()
rows = range(n)
for i in range(n):
    # This is i-th R^n basis vector
    ei.zero()
    ei[i] = 1
    # B.ei gives the i-th colum
    B.mult(ei, b)
    solver.solve(x, b)
    # Set the column
    col = [i]
    X.setValues(rows, col, x.array())
# Finalize
X.assemblyBegin()
X.assemblyEnd()
X = PETScMatrix(X)

# Check: (AX-B).y = 0
y = Vector(None, n)
y[:] = np.random.rand(n)
# Apply B to random y
By = Vector(None, n)
B.mult(y, By)
# AXy store to y
Xy = Vector(None, n)
X.mult(y, Xy)
A.mult(Xy, y)
# Error norm
y.axpy(-1, By)
print y.norm('linf')
