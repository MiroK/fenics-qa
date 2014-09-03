from dolfin import *
import numpy as np

nu = Constant(0.2)                  # poisson ratio (undamaged material)
u_f = 10                            # Total applied displacement
u_0 = 10**(-6)
ne = 100                            # Number of elements along one side of the domain
tol_b= 1E-14                        # Tolerence for defining the boundary

domain = Rectangle(-0.6, -0.2, 0.6, 0.2)
mesh = Mesh(domain, ne)

subdomains = CellFunction("size_t", mesh, 0)

class Omega0(SubDomain):
   def inside(self, x, on_boundary):
       return True if (-0.2 <= x[0] <= 0.2 and -0.1 <= x[1] <= 0.1) else False

subdomain0 = Omega0()
subdomain0.mark(subdomains, 1)

V0 = FunctionSpace(mesh, 'DG', 0)
E = Function(V0)

# Loop over all cell numbers, find corresponding subdomain number and fill cell value in E
E_values = [1, 0.1]                         # values of E in the two subdomains
for cell_no in range(len(subdomains.array())):
    subdomain_no = subdomains.array()[cell_no]
    E.vector()[cell_no] = E_values[subdomain_no]

plot(E)
interactive()
exit()

# V is a function space to solve the equilibrium equation over the domain
V = VectorFunctionSpace(mesh, 'Lagrange', 1)

 # Right boundary  (Applied displacement as tension)
def right_boundary(x, on_boundary):
    return on_boundary and abs(x[0]-0.6) < tol_b

def left_boundary(x, on_boundary):
    return on_boundary and abs(x[0]+0.6) < tol_b

# Strain Tensor
def eps(v):
    return sym(grad(v))

 #Definition of Stress Tensor components

def sigma11(v):
    return (eps(v)[0,0]+nu*eps(v)[1,1])

def sigma22(v):
    return (eps(v)[1,1]+nu*eps(v)[0,0])

def sigma12(v):
    return (1-nu)*eps(v)[0,1]
#-------------------------------------------------------------------------------------------------------
# bulk load (body forces)
f=Constant((0.0,0.0))
#-------------------------------------------------------------------------------------------------------
# Define the function, test and trial fields for equilibrium equation
u=TrialFunction(V)
u_n=TrialFunction(V)                   # Displacement field in current step
v =TestFunction(V)                     # Test function to solve equlibrium
# Bilinear form of elasticity
a= (E.vector()[cell_no]*(sigma11(u)*eps(v)[0,0]+ sigma22(u)*eps(v)[1,1]+sigma12(u)*eps(v)[0,1]))*dx
#a = (E.vector()[cell_no]*(eps(u)[0,0]*eps(v)[0,0]+nu*eps(u)[1,1]*eps(v)[0,0]+ nu*eps(u)[0,0]*eps(v)[1,1]+eps(u)[1,1]*eps(v)[1,1]+(1-nu)*eps(u)[0,1]*eps(v)[0,1]))*dx
# Linear form of elasticity
L=dot(f,v)*dx
T = 1                                # total simulation time
n_timesteps = 101                 # number of time steps
dt = T/float(n_timesteps-1)        # time step
#-------------------------------------------------------------------------------------------------------
times = np.linspace(0.,T,n_timesteps)
#for timestep in range(len(times)):
t = 0
while t <= T:# and Linf_dold < 1.0 :
    #plot(d_old, title="damage level")
    #t = times[timestep]linear
    t = t + dt
    print 'time =', t
    u_app = t*u_f
    print 'Applied disp', u_app
    u_n = Function(V)
    u_L = Constant((-1*u_app,0.0))
    u_R = Constant((u_app,0.0))
    Gamma_1 = DirichletBC(V, u_L, left_boundary)
    Gamma_2 = DirichletBC(V, u_R, right_boundary)
    bcs = [Gamma_1, Gamma_2]
    solve(a == L, u_n, bcs,solver_parameters={'linear_solver': 'cg','preconditioner': 'bjacobi'})
