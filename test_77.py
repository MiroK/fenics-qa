from dolfin import *
from math import sqrt
import numpy as np
import pylab

# Parameters
beta = 1.0/4.0                               # newmark-beta parameter
info("Newmark-beta = %f" % (beta))           
N = 100                                        # number of timesteps

rho = 1.6              # Density in g/cc 
mu = 10e-3          # Shear Modulus in GPA
inv = mu/rho      
c = sqrt(inv)          # Speed of wave (km/s)

#Define Mesh
dh  = 48
h = 1./dh                     # Mesh Size
mesh = UnitSquareMesh(dh,dh)

# Set time stepping information for CFL (mesh size - time step) condition
cfl_dt = h/c            
dt = 0.05*cfl_dt
info('cfl time step size = %.5f' % cfl_dt)
info('time step size = %.5f' % dt)
t = 0.0
T = N*dt

V = FunctionSpace(mesh,"CG",1)

def boundaries(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, Constant(0.0), boundaries)

ui = Expression('-40.*sin(pi*x[0])*sin(pi*x[1])*exp(-1000.*(pow(x[0]-0.5,2)+pow(x[1]-0.5,2)))')
u0 = interpolate(ui, V)                      # Initial conditions
p0 = interpolate(Constant(0.0),V)    #initial velocity

u = TrialFunction(V) 
v = TestFunction(V)  
u_new = Function(V)                   # displacement (solution)
p_new = Function(V)                   # velocity
a_new = Function(V)                   # acceleration

k = inner(nabla_grad(u),nabla_grad(v))*dx # stiffness matrix integrand
m = u*v*dx                                # mass matrix integrand

K = assemble(k) # assemble stiffness matrix
M = assemble(m) # assemble mass matrix

u0_vec = u0.vector()
p0_vec = p0.vector()

# Start the process by finding the initial acceleration M*a0 = f0 - K*u0
#a0.vector() = -K*u0_vec/M
#Define a function to solve the system or can I directly solve using the formula in the line above
def solveprob(m1, m2, u, bc): 
  prob = LinearVariationalProblem(m1, m2, u, bcs=bc)
  solver = LinearVariationalSolver(prob)
  solver.solve()

a0 = Function(V)                            #initial accceleration
m_2 = -u0*c**2*k
solveprob(m, m_2, a0, bc)

# Compute LU factorization of 'B' [[breaking down the terms from final discretised expression resulting from newmark beta]]
B = M + dt**2*K*beta*c**2        #LHS of discretised equation
bc.apply(B)
Bsolve = LUSolver(B)
Bsolve.parameters["reuse_factorization"] = True

while t <= T:

    # Perform predictor step 
    u_p = u0_vec + dt*p0_vec + dt*dt/2*(1-2*beta)*a0_vec    
    p_p = p0_vec + dt*a0_vec/2 
    L = -K*u_p  
    # Solve algebric system of equations
    bc.apply(L)
    Bsolve.solve(a_new.vector(), L)
    t += dt
    # Perform corrector step
    u_new.vector()[:] = u_p + dt**2*beta*a_new.vector()  
    p_new.vector()[:] = p_p + dt/2*a_new.vector()   
    u0.assign(u_new) # update old vector
    p0.assign(p_new)
