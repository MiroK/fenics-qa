from dolfin import *

A, B = 0.3, 0.5
middle = AutoSubDomain(lambda x: A - DOLFIN_EPS < x[0] < B + DOLFIN_EPS)

mesh = UnitIntervalMesh(1000)
cell_f = CellFunction('size_t', mesh, 0)
middle.mark(cell_f, 1)

V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

#dx = Measure('dx')[cell_f]
a = inner(Constant(1)*grad(u), grad(v))*dx
L = inner(Expression('sin(pi*x[0])'), v)*dx

bc_l = DirichletBC(V, Constant(1), 'x[0] < DOLFIN_EPS')
bc_r = DirichletBC(V, Constant(2), 'x[0] > 1- DOLFIN_EPS')
bc_c0 = DirichletBC(V, Constant(1.75), 'near(x[0], %g)' % A, 'pointwise')
bc_c1 = DirichletBC(V, Constant(1.75), 'near(x[0], %g)' % B, 'pointwise')
bcs = [bc_l, bc_r]#, bc_c0, bc_c1]

uh = Function(V)
solve(a == L, uh, bcs)

plot(uh)
interactive()

