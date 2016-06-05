from dolfin import *
mesh = UnitSquareMesh(200, 200)
V = FunctionSpace(mesh, "CG", 1)
u = TrialFunction(V)
v = TestFunction(V)
a = dot(grad(u), grad(v))*dx
m = u*v*dx

# Assemble the stiffness matrix and the mass matrix.
A = assemble(a)
M = assemble(m)
A = as_backend_type(A)
M = as_backend_type(M)

# Create eigensolver and compute eigenvalues
eigensolver = SLEPcEigenSolver(A)
eigensolver.parameters["spectrum"] = "smallest magnitude"
eigensolver.solve(3)

# Check the number of eigenvalues that converged.
# Extract the eigenvalues (ignore the imaginary part) and compare with the exact values
nconv = eigensolver.get_number_converged()
if MPI.rank(mesh.mpi_comm()) == 0:
    print "Number of eigenvalues successfully computed: ", nconv
eigenvalues = []
for i in range(nconv):
    r, c, rx, cx = eigensolver.get_eigenpair(i)  # real and complex part of eigenvalue
    eigenvalues.append(r)

if MPI.rank(mesh.mpi_comm()) == 0:
  print "Smallest positive eigenvalues computed and exact: "
  for i in range(min(nconv, 3)):
    print eigenvalues[i]
