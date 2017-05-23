from dolfin import *

mesh = UnitSquareMesh(40, 40)
f = Expression('sin(k*pi*x[1])+cos(k*pi*x[0])', k=5, degree=3)
g = Constant(1)
u0 = Constant(0)

V = FiniteElement('Lagrange', mesh.ufl_cell(), 1)
Q = FiniteElement('Real', mesh.ufl_cell(), 0)
W = MixedElement([V, Q])

W = FunctionSpace(mesh, W)
bc = DirichletBC(W.sub(0), u0, 'near(x[1], 0)')

u, p = TrialFunctions(W)
v, q = TestFunctions(W)

right = CompiledSubDomain('near(x[0], 1)')
bdries = FacetFunction('size_t', mesh, 0)
right.mark(bdries, 1)

ds1 = Measure('ds', domain=mesh, subdomain_data=bdries, subdomain_id=1)

a = inner(grad(u), grad(v))*dx + inner(p, v)*ds1 + inner(q, u)*ds1
L = inner(f, v)*dx + inner(g, q)*ds1

w = Function(W)
solve(a == L, w, bc)

u, _ = w.split(deepcopy=True)
plot(u)
interactive()

print assemble(u*ds1)
