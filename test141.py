from scipy.sparse import csc_matrix
from scipy.sparse.linalg import cg
from dolfin import *

parameters.linear_algebra_backend = "uBLAS"

class boundary(SubDomain):
   def inside(self,x,on_boundary):
      return on_boundary

uexact = Expression("x[0]*x[0] - x[1]*x[1]")

n = 100
mesh = UnitSquareMesh(n,n)
V = FunctionSpace(mesh, 'CG', 3)

v = TestFunction(V)
u = TrialFunction(V)

ubc = Function(V)

# Project boundary condition
M = assemble(u*v*ds)
b = assemble(uexact*v*ds)
bd = boundary()
bc = DirichletBC(V, ubc, bd)
binds = bc.get_boundary_values().keys()
rows, cols, values = M.data()
M = csc_matrix((values, cols, rows))
M = M[binds,:][:,binds]
rhs = b[binds]
x,info = cg(M, rhs)
ubc.vector().zero()
ubc.vector()[binds] = x

# Solve
a = inner(grad(u), grad(v))*dx
L = Constant(0)*v*dx
u = Function(V)
solve(a==L, u, bc)

# compute boundary stress du/dn
a = inner(grad(u), grad(v))*dx
b = assemble(a-L)
b = b[binds]
tau,info = cg(M, b)
stress = Function(V)
stress.vector().zero()
stress.vector()[binds] = tau

plot(u)
plot(stress)
interactive()
