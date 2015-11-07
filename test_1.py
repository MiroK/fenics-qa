from dolfin import *

# Create mesh and define function space
mesh = UnitSquareMesh(24, 24)
W1 = FunctionSpace(mesh, "CG", 1)
V = MixedFunctionSpace([W1,W1])


# Define test and trial functions
(u, ub) = TrialFunctions(V)
(v, w) = TestFunctions(V)


# Initial condition
c0 = Function(V)
s = Function(W1)
u_init = Expression("x[0]")
s.interpolate(u_init)
assign(c0.sub(0),s)
(u0, ub0) = as_vector((c0[0], c0[1]))

# Define parameters and coefficients
dt = 0.01
T = 0.1
k1 = 0.2
k2 = 0.1
D = 0.1

# Define variational problem
F1 = (1/dt)*(u-u0)*v *dx \
   + D*dot(grad(v), grad(u))*dx \
   - k1*u*v*dx + k2*ub*v*dx

F2 = (1/dt)*(ub-ub0)*v *dx \
    + D*dot(grad(w), grad(ub))*dx \
    - k2*ub*w*dx + k1*u*w*dx

F = F1 + F2
a = lhs(F); L = rhs(F)

# preassembly
A = assemble(a)
b = None

# Compute solution
c1 = Function(V)
t = dt

#plot(ub0); interactive()
file1 = File("solutions/u.pvd")
file2 = File("solutions/ub.pvd")
file3 = File("solutions/c.pvd")

# No need to create new function for each solution of the projection problem
comb0 = Function(W1)

# Function comb1 has the combination obtained by inteporpolating Expression
comb1 = Function(W1)

# Function comb3 will use assigner to get the combinations
comb2 = Function(W1)
# Assigner from first component to W1
first_assigner = FunctionAssigner(W1, V.sub(0))
# Assigner from second component to W1
second_assigner = FunctionAssigner(W1, V.sub(1))
comb_aux = Function(W1)  # Auxiliary function


while t < T:
    #Write to file
    (u0, ub0) = c0.split()
    file1 << u0
    file2 << ub0

    #Solve problem
    b = assemble(L, tensor=b)
    solve(A, c1.vector(), b)
    c0.assign(c1)
    (u1, ub1) = c1.split()

    # Sum solutions: c = u + ub
    c = TrialFunction(W1)
    rho = TestFunction(W1)
    L1 = ub1*rho*dx + u1*rho*dx
    a1 = c*rho*dx
    # Solve into comb0
    solve(a1 == L1, comb0)
    file3 << comb0
    
    # Combine u1 and ub1 via interpolation of Expression
    comb1.interpolate(Expression('f1+f2', f1=u1, f2=ub1))
    # Check the difference
    comb1.vector().axpy(-1, comb0.vector())
    print('\t comb1 error %g' % comb1.vector().norm('l2'))
    
    first_assigner.assign(comb_aux, u1)   # comb_aux has u1
    second_assigner.assign(comb2, ub1)   # comb2 has ub1 
    comb2.vector().axpy(1, comb_aux.vector())   # add comb_aux to comb2
    # Check the difference
    comb2.vector().axpy(-1, comb0.vector())
    print('\t comb2 error %g' % comb2.vector().norm('l2'))
    

    # print("Solution vector norm (0): {!r}".format(combined.vector().norm("l2")))
    t += dt


# Save last solution to file
file1 << u0
file2 << ub0
