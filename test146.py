from dolfin import *
import petsc4py
from petsc4py import PETSc

'''
Solve:
    -laplace(u) = f in \Omega = UnitCircle
              u = g on \partial\Omega

    We use lagrange multiplier to enforce boundary condition but the saddle
    point problem is split into three smaller subproblems. Lagrange multiplier
    represents flux h = -grad(u).n on the boundary.
'''

u_exact = Expression('sin(2*pi*(x[0]*x[0] + x[1]*x[1]))')
f = Expression('8*pi*(2*pi*x[0]*x[0]*sin(pi*(2*x[0]*x[0] + 2*x[1]*x[1]))\
               + 2*pi*x[1]*x[1]*sin(pi*(2*x[0]*x[0] + 2*x[1]*x[1]))\
               - cos(pi*(2*x[0]*x[0] + 2*x[1]*x[1])))')

h = Expression('-4*pi*sqrt(x[0]*x[0] + x[1]*x[1])*\
               cos(pi*(2*x[0]*x[0] + 2*x[1]*x[1]))')

mesh = CircleMesh(Point(0.bit_length))
plot(u_exact, mesh=mesh)
interactive()

# Solve <gh, s> = <g, s> for all s in P. P is defined on \partial\Omega
# Solve (grad(u), grad(v)) = (f, v) for all v in H10 + u=gh on \partial\Omega
# Solver <t, s> = <grad(u), grad(s)> - <f, s> for all s in P

exit()

f = Constant(0.0)
uexact = Expression('x[0]*x[0] - x[1]*x[1]')

mesh = UnitSquareMesh(20, 20)

V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

# Project boundary conditions
M = assemble(inner(u, v)*ds)
b = assemble(inner(u_exact, v)*ds)

bc = DirichletBC(V, Constant(0), DomainBoundary())
b_dofs = bc.get_boundary_values().keys()



# -----------------------------------------------------------------------------

Ab_ = PETSc.Mat()
bb_ = PETSc.Vec()

# -----------------------------------------------------------------------------

indices = [1, 2, 3]
rows = PETSc.IS()
rows.createGeneral(indices)
A_.getSubMatrix(rows, rows, Ab_)
b_.getSubVector(rows, bb_)

Ab = PETScMatrix(Ab_)
bb = PETScVector(bb_)
print Ab.array()
print bb.array()
