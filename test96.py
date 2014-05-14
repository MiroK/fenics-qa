from dolfin import *
import math
import ufl

mesh = UnitCubeMesh(1,1,1)

V = VectorFunctionSpace(mesh, 'Lagrange', 2) # displacement
Q = FunctionSpace(mesh,'Lagrange',1) # pressure
R = FunctionSpace(mesh,'Lagrange',2) # bulk concentration
M = FunctionSpace(mesh,'Lagrange',2) # Lagrange multiplier

W = MixedFunctionSpace([V,Q,R,M])

# Test and trial functions
vqrm  = TestFunction(W) # Test function in the mixed space
v,q,r,tcs = split(vqrm)
dupccs = TrialFunction(W) # Trial function in the mixed space (solution increment)

upccs = Function(W) # Function in the mixed space to store the solution
u,p,c,cs = split(upccs) # Function in each subspace to write the functional

# nonlinear effective diffusion coefficient
def g(c):
    return D*((2*chi-1.)*Omega*Jo*c-1.)/pow((1.+Omega*Jo*c),3)

def dmudc(c):
    return -(Rg*T)*((2.*chi-1.)*Omega*Jo-1./c)/pow((1.+Omega*Jo*c),3)

# Parameters
Omega = 6.023e-5
D = 1E-1
Gdry = 4E5
Rg = 8.314
T = 293.0
lamo = 1.001
Jo = lamo**3
cod = (lamo**3-1.0)/Omega
co = cod/Jo
chi = 0.2
po = Gdry/lamo
muo = Rg*T*(math.log(Omega*cod/(1+Omega*cod))+1/(1+Omega*cod)+ chi/pow((1+Omega*cod),2)) \
+Omega*po
mus = muo
mutop = -10000

# Kinematics
I = Identity(3)
F = I + grad(u)
J = det(F)

# Constitutive equations

S = (Gdry/lamo)*F-p*cofac(F)
mu = Rg*T*(ufl.ln(Omega*Jo*c/(1+Omega*Jo*c)) + 1./(1.+Omega*Jo*c)+ \
         chi/pow((1.+Omega*Jo*c),2))+Omega*p

h = g(c)*grad(c)+(-c*D/(Rg*T))*Omega*grad(p)

dmudp = Omega
dmu = mu-mus

# Boundary conditions

def left_boundary(x, on_boundary): # x = 0
    return on_boundary and abs(x[0]) < DOLFIN_EPS

def back_boundary(x, on_boundary): # y = 0
    return on_boundary and abs(x[1]) < DOLFIN_EPS

def bottom_boundary(x, on_boundary): # z = 0
    return on_boundary and abs(x[2]) < DOLFIN_EPS

bcl = DirichletBC((W.sub(0)).sub(0), Constant(0.), left_boundary) # u = 0 on x = 0
bcb = DirichletBC((W.sub(0)).sub(1), Constant(0.), back_boundary) # v = 0 on y = 0
bcbo = DirichletBC((W.sub(0)).sub(2), Constant(0.), bottom_boundary) # w = 0 on z = 0
bc = [bcl,bcb,bcbo]

# WEAK FORM OF BALANCE EQUATIONS

a = inner(S,grad(v))*dx + q*(J-(1./Jo+Omega*c))*dx + inner(h,grad(r))*dx + inner(dmu,tcs)*ds \
        + inner(dmudc(c)*cs,r)*ds + inner(dmudp*cs,q)*ds
Jac = derivative(a,upccs,dupccs)

problem = NonlinearVariationalProblem(a, upccs, bc, Jac)
solver = NonlinearVariationalSolver(problem)

# Initial guess for unknowns
c = interpolate(Constant(co), R)
p = interpolate(Constant(po), Q)
cs = interpolate(Constant(co),M)

# solve non-linear problem
solver.solve()

# plot solvent concentration
w = Function(V)
w = upccs.split()[2]
plot(w)
interactive()
