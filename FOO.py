from petsc4py import PETSc
from dolfin import *

uexact = Expression('x[0]*x[0] - x[1]*x[1]')

n = 40
mesh = UnitSquareMesh(n, n)
V = FunctionSpace(mesh, 'CG', 1)

u = TrialFunction(V)
v = TestFunction(V)
ubc = Function(V)

# Project boundary condition
M = PETScMatrix()
b = PETScVector()
assemble(u*v*ds, tensor=M)
assemble(uexact*v*ds, tensor=b)
bc = DirichletBC(V, ubc, DomainBoundary())
binds = bc.get_boundary_values().keys()
binds.sort()

# Extract submatrices
bdr_dofs = PETSc.IS()
bdr_dofs.createGeneral(binds)
M_petsc = PETSc.Mat()
b_petsc = PETSc.Vec()
x_petsc = PETSc.Vec()
M.mat().getSubMatrix(bdr_dofs, bdr_dofs, M_petsc)
b.vec().getSubVector(bdr_dofs, b_petsc)
b.vec().getSubVector(bdr_dofs, x_petsc)
M = PETScMatrix(M_petsc)
rhs = PETScVector(b_petsc)
x = PETScVector(x_petsc)
x.zero()
print M.size(0), M.size(1)

# Solve the projection system
solve(M, x, rhs)

# Assign x to boundary vector
ubc.vector()[binds] = x

# Solve
a = inner(grad(u), grad(v))*dx
L = Constant(0)*v*dx
u = Function(V)
solve(a==L, u, bc)

# compute boundary stress du/dn
a = inner(grad(u), grad(v))*dx
assemble(a-L, tensor=b)
b.vec().getSubVector(bdr_dofs, b_petsc)
b = PETScVector(b_petsc)
solve(M, x, b)
stress = Function(V)
stress.vector()[binds] = x

plot(u)
plot(stress)
interactive()
