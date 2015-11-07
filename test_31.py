from dolfin import *

left = AutoSubDomain(lambda x: x[0] < 0 + DOLFIN_EPS)
right = AutoSubDomain(lambda x: x[0] > 0 - DOLFIN_EPS)

f_left = Constant(1.)
f_right = Constant(2.)

kappa_left = Constant(1.)
kappa_right = Constant(10.)

mesh = IntervalMesh(100, -1, 1)
cell_f = CellFunction('size_t', mesh, 0)
left.mark(cell_f, 1)
right.mark(cell_f, 2)

dx_left = dx(1, domain=mesh, subdomain_data=cell_f)
dx_right = dx(2, domain=mesh, subdomain_data=cell_f)

V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)
bc = DirichletBC(V, Constant(0), DomainBoundary())

a = inner(kappa_left*grad(u), grad(v))*dx_left + \
    inner(kappa_right*grad(u), grad(v))*dx_right

L = inner(f_left, v)*dx_left + inner(f_right, v)*dx_right

uh = Function(V)
solve(a == L, uh, bc)

plot(uh)
interactive()

