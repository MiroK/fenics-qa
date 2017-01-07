from fenics import *
import numpy as np
import matplotlib.pyplot as plt

# Parameters
c = 1.25
l = 5.0 
P_surf = 0.0
sigma_0 = 3.0 
gamma = 0.65 
P0 = sigma_0*gamma 
Kv = 10 
cm = 5
# Derivative of pressure at the bottom surface.
der_p_bottom = 0.0
# Derivative of displacement at the top surface.
der_w_top = sigma_0/Kv

# The initial time.
t = 0.0
# The final time
T = 2.0
# The number of time steps to take.
num_steps = 40
# The step size
dt = (T-t)/num_steps


# Create mesh and define function space.
mesh = IntervalMesh(100,0,l)
P1 = FiniteElement('P', interval, 1)
element = MixedElement([P1, P1])
V = FunctionSpace(mesh, element)

# Define boundary conditions for pressure.
tol = 1E-14

def upper_boundary(x,on_boundary):
    return on_boundary and abs(x[0])<tol

# An expression for the derivative of P at the bottom surface.
der_p_D = Expression('der_p_bottom*x[0]/l',degree = 1, \
    der_p_bottom=der_p_bottom, l=l)

bc_p = DirichletBC(V.sub(0),Constant(P_surf),upper_boundary)

# Define boundary conditions for displacements.
def bottom_boundary(x,on_boundary):
    return on_boundary and abs(x[0]-l)<tol

der_w_D = Constant(der_w_top)

bc_w = DirichletBC(V.sub(1),Constant(0.0),bottom_boundary)

# Aggregate the Dirichlet boundary conditions in a list.
bc = [bc_p, bc_w]

# Define the test function space.
v_1, v_2 = TestFunctions(V)

# Define functions for solutions at previous and current time steps.
# Previous time step.
# Define the initial conditions.
init_cond = Expression(('P0','sigma_0/Kv*(l-x[0])'), degree = 1,
                       sigma_0=sigma_0, Kv=Kv, l=l,P0=P0)
u_n = interpolate(init_cond,V)
p_n,w_n = u_n.split()            
# Current time step.
u = Function(V)
p_, w_ = split(u)

# Define the variational problem.
F = p_*v_1*dx + c*dt*inner(grad(p_),grad(v_1))*dx - c*dt*der_p_D*v_1*ds \
    - p_n*v_1*dx \
    + inner(grad(w_),grad(v_2))*dx + cm*inner(grad(p_),grad(v_2))*dx-der_w_D*v_2*ds

# Change the variational problem into a bilinear form.
a,L = lhs(F),rhs(F)

# Now, declare the current time step's variable as a function.
p_ = Function(V)

# Create the progress bar.
progress = Progress('Time-stepping')
set_log_level(PROGRESS)

# Open the file to store results in VTK format.
vtkfile_w = File('Terzaghi_results_FI/pressure.pvd')

# Form and solve the linear system.
for n in range(num_steps):

    # Update current time.
    t += dt

    # Solve the variational problem.
    solve(a == L, u, bc)

    # Extract the solution.
    p_,w_ = u.split()
    p_comp = p_.compute_vertex_values(mesh)
    w_comp = w_.compute_vertex_values(mesh)

    # Write the results to a VTK file for visualization in Paraview.
    vtkfile_w << (p_,t)
    vtkfile_w << (w_,t)

    # Update the previous solution.
    u_n.assign(u)

    # Update the progress bar.
    progress.update(t/T)
