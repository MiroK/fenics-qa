from dolfin import *

# Create mesh and space
mesh = UnitIntervalMesh(2)
V = FunctionSpace(mesh, 'CG', 1)
dim = V.dim()

# Setup zero matrix
Mat = PETScMatrix()
mat = Mat.mat()
mat.create()
mat.setSizes((dim, dim))
mat.setType('aij')
mat.setUp()
# Add some entries
for p in range(dim):
     mat[p, p] = 1.
mat.assemble()
A1 = PETScMatrix(mat)
print A1.array()

# Define second matrix
u = TrialFunction(V)
v = TestFunction(V)
a = inner(grad(u), grad(v))*dx
A2 = PETScMatrix()
assemble(a, tensor=A2)
print A2.array()

# Add matrices in order to solve the associated linear system
A2.axpy(1, A1, False)
print A2.array()
# A3 = A1 + as_backend_type(A2)
