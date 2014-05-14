from dolfin import *

NUM_CELL = 24
mesh = UnitSquareMesh(NUM_CELL,NUM_CELL)
h = CellSize(mesh)
n = FacetNormal(mesh)

# Create FunctionSpaces
V =  FunctionSpace(mesh, "CG", 3)
W =  FunctionSpace(mesh, "CG", 1)
DG = FunctionSpace(mesh, "DG", 0)

# Boundary conditions
def right(x, on_boundary): return x[0] > (1.0 - DOLFIN_EPS)
def left(x, on_boundary): return x[0] < DOLFIN_EPS
def bottom_center(x, on_boundary):
    return x[1] < DOLFIN_EPS and (x[0] > 1./3. - DOLFIN_EPS and x[0] < 2./3.0 + DOLFIN_EPS)
def bottom_lr(x, on_boundary):
    return x[1] < DOLFIN_EPS and (x[0] < 1./3. + DOLFIN_EPS or x[0] > 2./3.0 - DOLFIN_EPS)
def top(x, on_boundary):
    return x[1] > 1.0 - DOLFIN_EPS

g0 = Constant(0.0)
g1 = Constant(1.0)

bc0 = DirichletBC(V, g1, bottom_center)
bc1 = DirichletBC(V, g0, bottom_lr)
bc3 = DirichletBC(V, g0, top)
bc4 = DirichletBC(V, g0, right)

bcs = [bc0, bc1, bc4, bc3]

# Parameters
epsilon = Constant(0.000000001)
c = Constant(0.)
b = Expression(('-x[1]', 'x[0]'))
f = Constant(0.)

uh = TrialFunction(V)
vh = TestFunction(V)
wh = TestFunction(W)
dg = TestFunction(DG)

# SUPG (SDFEM) method
bb = assemble(dot(b,b)*wh*dx)
bF = Function(W,bb)
tau = 1.
a1 = (epsilon*dot(grad(uh),grad(vh)) + vh*dot(b,grad(uh)) + c*uh*vh)*dx
a2 = (h/(2.*sqrt(bF))*tau*inner(dot(b,grad(uh)),dot(b,grad(vh))))*dx
a = a1 + a2
L = f*vh*dx + h/(2.*sqrt(bF))*tau*h*f*dot(b, grad(vh))*dx

# Compute solution
uh = Function(V)
solve(a == L, uh, bcs)

# Error estimators
residual = dg*(-epsilon*div(grad(uh))+dot(b,grad(uh))+c*uh-f)**2*dx
indicators = assemble(residual)


# HOW TO PLOT INDICATORS?
print indicators.min(), indicators.max()
plot(Function(DG, indicators)) # here it fails



error_estimate = sqrt(sum(i for i in indicators.array()))
print "error_estimate = ", error_estimate

# Plot solution and mesh
plot(uh)

# Hold plot
interactive()
