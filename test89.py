from dolfin import *

# Create mesh and define function space
m_num = 6 # number of mesh in each dim
mesh = UnitCubeMesh(m_num, m_num, m_num)
dim = mesh.topology().dim() 

##  Function Spaces
Pu = VectorFunctionSpace(mesh, "Lagrange", 1) # space for displacements
Pp = FunctionSpace(mesh, "Lagrange", 1)     # space for pressure
V = MixedFunctionSpace([Pu,Pp])                    # mixed space

# Functions
dup = TrialFunction(V)    # Incremental displacement-pressure
wq  = TestFunction(V)     # Test function
up  = Function(V)         # Displacement-pressure from previous iteration
u, p = up.split()          # Function in each subspace to write the functional
w, q = wq.split()         # Test Function split

# Loads
B  = Constant((0.0,  0.0, 0.0))  # Body force per unit volume
Trac = Expression(("0.0","0.0","1.0")) # Traction force on the boundary
g_bar = Constant(0.0)            # Normal flux


## Boundary Condition

# Create mesh function over cell facets
exterior_facet_domains = FacetFunction("size_t", mesh)

# Mark Neumann boundaries
class TopBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14   # tolerance for coordinate comparisons
        return on_boundary and abs(x[2] - 1.0) < tol

Gamma_T = TopBoundary()
exterior_facet_domains.set_all(1)
Gamma_T.mark(exterior_facet_domains, 0)
ds_neumann = ds[exterior_facet_domains]

# Mark Dirichlet boundaries
top    = CompiledSubDomain("near(x[2], side) && on_boundary", side = 1.0)
bottom = CompiledSubDomain("near(x[2], side) && on_boundary", side = 0.0)

# Assign Dirichlet boundaries (x = 0 or x = 1)
d_top    = Expression(("0.0", "0.0", "0.0"))
d_bottom = Expression(("0.0", "0.0", "0.0"))
bc_top    = DirichletBC(V.sub(0), d_top, top)
bc_bottom = DirichletBC(V.sub(0), d_bottom, bottom)

p_top    = Expression("0.0")
p_bottom = Expression("1.0")
bc_ptop = DirichletBC(V.sub(1), p_top, top)
bc_pbottom = DirichletBC(V.sub(1), p_bottom, bottom)

bcs = [bc_bottom, bc_ptop]

## Initial conditions
up_1 = Function(V)
u_1.interpolate(Constant((0.0, 0.0, 0.0, 0.0, 0.0)))
u_1, p_1 = up1_1.split()

## Kinematics
I = Identity(dim)             # Identity tensor
F = I + grad(u)             # Deformation gradient - seems like underyling coordinates is initial
F = variable(F)             # Make F a variable for tensor differentiations
C = F.T*F                   # Right Cauchy-Green tensor
invF = inv(F)

# virtual Kinematics
ddotF = grad(w)
ddotE = 1.0/2.0 * (ddotF.T*F + F.T*ddotF)

# Invariants of deformation tensors
J  = det(F)
Ic = tr(C)
IIIc = det(C)

# Elasticity parameters
Ee, nu = 10.0, 0.45
mu, lmbda = Constant(Ee/(2*(1 + nu))), Constant(Ee*nu/((1 + nu)*(1 - 2*nu)))

# Permeability
perm = 0.1
K_perm = perm*I

## Potential Energy
# Strain energy density (compressible neo-Hookean model)
psi = (mu/2.0)*(Ic - 3) - mu*ln(J) + (lmbda/2.0)*(ln(J))**2

# PK1 stress tensor
P = diff(psi,F)
# PK2 stress tensor
S = inv(F)*P

## time steps
dt = 0.1 # time step
T_toal = 1.0
t = dt
tn=0

# Compute residual
R = (inner(S, ddotE) + dot(J*(K_perm*invF.T*grad(p)), invF.T*grad(q)))*dx \
    - p*J*inner(ddotF, invF.T)*dx \
    + (q*J*inner(grad((u-u_1)/dt),invF.T))*dx \
    + (inner(B,w))*dx - (inner(Trac,w))*ds_neumann(0) \
    + (inner(g_bar,q))*ds_neumann(1)

# Compute Jacobian of F
Jac = derivative(R, up, dup)

# Set up the problem
problem = NonlinearVariationalProblem(R, up, bcs=bcs, J=Jac )
solver = NonlinearVariationalSolver(problem)
solver.parameters["newton_solver"]["linear_solver"] = "bicgstab"
solver.parameters["newton_solver"]["preconditioner"] = "ilu"

## Run Simulation
while t<T_toal:
    print 'time = ', t    

    # solve
    solver.solve()

    # update
    t += dt
    tn+=1
    up_1.assign(up)

    plot(u_1, interactive=True)
    plot(p_1, interactive=True)
