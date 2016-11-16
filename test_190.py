from dolfin import *

mesh = UnitCubeMesh(5, 5, 5)
V = VectorFunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

kappa, mu = Constant(2), Constant(1)

eps = variable(sym(grad(u)))
psi =  kappa / 2. * tr(eps) ** 2 + mu * inner(dev(eps), dev(eps))
sigma = diff(psi, eps)

# Uncommenting this yields the correct solution
eps1 = sym(grad(u))
sigma1 = kappa * tr(eps1) * Identity(len(u)) + 2.0 * mu * dev(eps1)

a = inner(sigma, sym(grad(v)))*dx
a1 = inner(sigma1, sym(grad(v)))*dx

A, A1 = map(lambda form: as_backend_type(assemble(form)).mat(), (a, a1))
A.axpy(-1., A1)
print A.norm(0)
