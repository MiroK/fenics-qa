from dolfin import *

mesh = UnitSquareMesh(10, 10)

subdomains = CellFunction('size_t', mesh, 0)
for cell in cells(mesh):
    M = cell.midpoint()
    if M.y() > M.x():
        subdomains[cell] = 1

plot(subdomains)

V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

f0 = Constant(2)  # source to be used in subdomain 0
f1 = Constant(1)  # source to be usde in subdomain 1

dx = Measure('dx')[subdomains]
a = inner(f1*grad(u), grad(v))*dx(0) + inner(f0*grad(u), grad(v))*dx(1)
# Inluclude sources
L = inner(f0, v)*dx(0) + inner(f1, v)*dx(1)
bc = DirichletBC(V, Constant(0), DomainBoundary())
uh = Function(V)

solve(a == L, uh, bc)
plot(uh)
interactive()
