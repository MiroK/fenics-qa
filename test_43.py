from dolfin import *

r0 = 1
r1 = 2
mesh = IntervalMesh(100, r0, r1)

V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

bc1 = DirichletBC(V, Constant(0), near(x[0], r1))
bc0 = DirichletBC(V, Constant(1), near(x[0], r0))
bcs = [bc0, bc1]

a = inner(grad(u), grad(v))*Expression('x[0]')*dx


