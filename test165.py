from dolfin import *

# Constants
ys = 1
t = 0
dt = 0.01
T = 10
rho1 = 1
rho2 = 10
mu1 = 1
mu2 = 0
lamda1 = 1
lamda2 = 1

# Mesh, functionspace, functions
mesh = RectangleMesh(0,0,2,2,20,20)
V = VectorFunctionSpace(mesh, "CG", 1)
D = FunctionSpace(mesh, "DG", 0)
u = TrialFunction(V)
v = TestFunction(V)

# Define Subdomains
cf = CellFunction("size_t", mesh, 0)
domain1 = AutoSubDomain(lambda x: x[1] > ys)
domain1.mark(cf,1)

# Define variable constants
rhof = Expression("x[1] > ys ? rho1 : rho2",
                  ys=ys,rho1=rho1,rho2=rho2)
muf = Expression("x[1] > ys ? mu1 : mu2",
                 ys=ys,mu1=mu1,mu2=mu2)
lamdaf = Expression("x[1] > ys ? lamda1 : lamda2",
                    ys=ys,lamda1=lamda1,lamda2=lamda2)
stepf = Expression("x[1] > ys ? pow(dt,2)/rho1 : pow(dt,2)/rho2",
                   ys=ys,dt=dt,rho1=rho1,rho2=rho2)

rho = interpolate(rhof, D)
mu = interpolate(muf, D)
lamda = interpolate(lamdaf, D)
step = interpolate(stepf, D)

# Stress tensor
def sigma(u, lamda, mu):
    return lamda*div(u)*Identity(2) + mu*(grad(u) + grad(u).T)

# Initial conditions
Ixy = Constant((1,0))
Vxy = Constant((2,0))
u2 = interpolate(Ixy, V)
u1 = interpolate(Vxy, V)

# Variational forms
F = inner(u,v)*dx - 2*inner(u1,v)*dx + inner(u2,v)*dx + \
    step*inner(sigma(u1,lamda,mu),grad(v))*dx

t = 2*dt
left = assemble(lhs(F))
u = Function(V)

# Time stepping
while t < T + DOLFIN_EPS:
    plot(u2, rescale=False)
    right = assemble(rhs(F))
    begin("solving at time step t=%g" % t)
    solve(left, u.vector(), right)
    end()

    u2.assign(u1)
    u1.assign(u)
