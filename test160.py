from dolfin import *

# Mesh
mesh = UnitSquareMesh(10,10)

# FunctionSpaces on the mesh
V =   FunctionSpace(mesh, "CG", 1)
DG0 = FunctionSpace(mesh, "DG", 0)

# Boundary conditions
def whole_boundary(x, on_boundary):
  return on_boundary
bcs = DirichletBC(V, 0., whole_boundary)

# Data
epsilon = Constant(1.e-8)
c = Constant(0.)
b = Constant((1., 0.))
b_perp = Constant((0., 1.)) # define b_perp automatically based on b???
f = Constant(1.)
tau = Constant(0.1)

# Trial And Test Functions
u = TrialFunction(V)
v = TestFunction(V)
dg0 = TestFunction(DG0)

# Solve for uh using SDFEM
a = (epsilon*dot(grad(u),grad(v)) + v*dot(b,grad(u)) + c*u*v)*dx +\
    inner(-epsilon*div(grad(u))+dot(b,grad(u))+c*u,tau*dot(b,grad(v)))*dx
L = f*v*dx + inner(f,tau*dot(b,grad(v)))*dx
uh = Function(V)
solve(a == L, uh, bcs)

# Compute a value of an indicatodef fcn_in_ind(x):def fcn_in_ind(x):
fcn_in_ind = lambda u: conditional(gt(u, 1), sqrt(u), 2.5*u**2 - 1.5*u**3)

indicator_ufl = ((-epsilon*div(grad(uh))+dot(b,grad(uh))+c*uh-f)**2 +
  fcn_in_ind(abs(dot(b_perp,grad(uh)))))*dg0*dx # use fcn_in_ind() in the ufl form instead of sqrt(...)???
indicator_assemble = assemble(indicator_ufl)
error_estimate = sum(i for i in indicator_assemble)
print "error_estimate = ", error_estimate

# Plot solution and mesh
plot(uh)

# Hold plot
interactive()


