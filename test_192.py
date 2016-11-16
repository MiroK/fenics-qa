from __future__ import print_function
from dolfin import *

# Define mesh, function space
mesh = UnitSquareMesh(32, 32)
V = FunctionSpace(mesh, "Lagrange", 1)

# Define basis and bilinear form
u = TrialFunction(V)
v = TestFunction(V)
a = dot(grad(u), grad(v))*dx

# Assemble stiffness form
A = PETScMatrix()
#The below statement is OK.
#assemble(a, tensor=A)
#But using the below command replacing assemble(a,tensor=A) will cause MPI_Comm faulty.
A=as_backend_type(assemble(a))
print(type(A))

# Create eigensolver
eigensolver = SLEPcEigenSolver(A)

# Compute all eigenvalues of A x = \lambda x
eigensolver.solve()

# Extract largest (first) eigenpair
r, c, rx, cx = eigensolver.get_eigenpair(0)
print ("Largest eigenvalue: ", r,c)
