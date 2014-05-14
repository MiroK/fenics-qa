from dolfin import *

mesh = UnitSquare(12,12)

# Create FunctionSpaces
V =  FunctionSpace(mesh, 'CG', 1)

# Boundaries case
def dOmega(x, on_boundary): return x[0] > (1.0 - DOLFIN_EPS) or x[0] < DOLFIN_EPS or x[1] > (1.0 - DOLFIN_EPS) or x[1] < DOLFIN_EPS
g0 = Constant(0.0)

# No-slip boundary condition for velocity
bc0 = DirichletBC(V, g0, dOmega)

# Parameters and functions
f = Constant(1.)
uh = TrialFunction(V)
vh = TestFunction(V)
a = dot(grad(uh), grad(vh))*dx
L = f*vh*dx

# Compute solution
uh = Function(V)
solve(a == L, uh, bc0)

# Mark boundary adjacent cells
boundary_adjacent_cells = [myCell for myCell in cells(mesh)
                                  if any([myFacet.exterior() for myFacet in facets(myCell)])]
cell_domains = CellFunction('size_t', mesh)
cell_domains.set_all(1)
for myCell in boundary_adjacent_cells:
    cell_domains[myCell] = 0

# Plot cell_domains
plot(cell_domains, interactive=True, title='cd')

 
dx = Measure('dx')[cell_domains]
integral = assemble(uh*(dx(0)+dx(1)))
# Try to integrate over the whole domain
print 'Integral = %e\r'%integral
# and over interior cells
integral2 = assemble(uh*dx(1))
print 'Integral2 = %e\r'%integral2

# Plot solution and mesh
plot(uh)
plot(mesh)

# Hold plot
interactive()
