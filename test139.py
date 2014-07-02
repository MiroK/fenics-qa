from dolfin import *

# Create mesh and define function space
#geom = Rectangle(0., 0., 5., 5.) - Rectangle(0., 0., 2, 2)
#diff0 = Rectangle(0, 0, 3, 3)
#diff1 = Rectangle(0:, 0, 3, 3) - Rectangle(0, 0, 2, 2)
#geom.set_subdomain(1, diff1)

#geom = Rectangle(0., 0., 5., 5.) - Rectangle(0., 0., 2, 2)
#geom.set_subdomain(1, Rectangle(0, 0, 3, 3))

geom = Rectangle(0., 0., 5., 5.) - Circle(0., 0., 2.5)
diff = CSGIntersection(Circle(0., 0., 3.5) - Circle(0., 0., 2.51),
                       Rectangle(0.01, 0.01, 3.5, 3.5))
geom.set_subdomain(1, diff)
mesh = Mesh(geom, 35)

plot(mesh)
V = FunctionSpace(mesh, 'Lagrange', 1)

print MeshQuality.radius_ratio_min_max(mesh)

# Define boundary conditions
u0 = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]')

#def u0_boundary(x, on_boundary):
#    return on_boundary
f_f = FacetFunction('size_t', mesh, 0)
DomainBoundary().mark(f_f, 1)
plot(f_f)
plot(u0, mesh=mesh)
bc = DirichletBC(V, u0, DomainBoundary())

interactive()
# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = inner(nabla_grad(u), nabla_grad(v))*dx('everywhere')
L = f*v*dx('everywhere')

# Compute solution
u = Function(V)
A = assemble(a)
b = assemble(L)

from collections import defaultdict
import numpy as np

dofs = []
xs = V.dofmap().tabulate_all_coordinates(mesh).reshape((-1, 2))
for x in xs:
    if np.hypot(x[0], x[1]) < DOLFIN_EPS:
        dofs.append('0')
    else:
        dofs.append(' '.join([str(x[0]), str(x[1])]))

print 'Total dofs =', len(dofs), 'Unique dofs =', len(set(dofs))
print xs

bc.apply(A, b)
solve(A, u.vector(), b)

# Plot solution and mesh
plot(u)

# Dump solution to file in VTK format
file = File('poisson.pvd')
file << u

# Hold plot
