from dolfin import *

# Load mesh and subdomains
mesh = Mesh("dolfin_fine.xml.gz")
sub_domains = MeshFunction("size_t", mesh, "dolfin_fine_subdomains.xml.gz")

# Define function spaces
V1 = VectorFunctionSpace(mesh, 'CG', 1)
V2 = VectorFunctionSpace(mesh, 'Bubble' ,3)
Q = FunctionSpace(mesh, 'CG', 1)

u, b, p = map(TrialFunction, (V1, V2, Q))
v, c, q = map(TestFunction, (V1, V2, Q))

# No-slip boundary condition for velocity
# x1 = 0, x1 = 1 and around the dolphin
noslip = Constant((0, 0))
bc0 = DirichletBC(V1, noslip, sub_domains, 0)

# Inflow boundary condition for velocity
# x0 = 1
inflow = Expression(("-sin(x[1]*pi)", "0.0"), degree=2)
bc1 = DirichletBC(V1, inflow, sub_domains, 1)

# Upper triangular part of A
a00 = inner(grad(u), grad(v))*dx
a01 = inner(grad(b), grad(v))*dx
a02 = inner(p, div(v))*dx

a10 = inner(grad(c), grad(u))*dx
a11 = inner(grad(b), grad(c))*dx
a12 = inner(p, div(c))*dx

a20 = inner(q, div(u))*dx
a21 = inner(q, div(b))*dx 

f = Constant((0, 0))
# Rhs components
L0 = inner(f, v)*dx
L1 = inner(f, c)*dx

from block import block_assemble, block_bc, block_mat
from block.iterative import MinRes
from block.algebraic.petsc import ML, Cholesky

# lhs
AA = block_assemble([[a00, a01, a02],
                     [a10, a11, a12],
                     [a20, a21,   0]])
# rhs
b = block_assemble([L0, L1, 0])

# Collect boundary conditions
bcs = [[bc0, bc1], [], []]
block_bc(bcs, True).apply(AA).apply(b)

b22 = inner(p, q)*dx
B22 = assemble(b22)
# Possible action of preconditioner
BB = block_mat([[ML(AA[0, 0]), 0, 0],
                [0, Cholesky(AA[1, 1]), 0],
                [0, 0, Cholesky(B22)]])

AAinv = MinRes(AA, precond=BB, tolerance=1E-10, maxiter=500)
U, B, P = AAinv*b

p = Function(Q)
p.vector()[:] = P

# Plot solution
plot(p)
interactive()
