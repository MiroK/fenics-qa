from mshr import Sphere, generate_mesh
from dolfin import *
import numpy as np

# Radius of outer and inner sphere
oradius, iradius = 5., 1.
tol = 1e-3

#Material parameters
nu = 0.3
mu = 1.0
Young = 2.*mu*(1.+nu)
lmbda = 2.*mu*nu/(1.-2.*nu)

# Geometry
outer_sphere = Sphere(Point(0., 0., 0.), oradius)
inner_sphere = Sphere(Point(0., 0., 0.), iradius)
g3d = outer_sphere - inner_sphere
mesh = generate_mesh(g3d, 8)

#Define the outer sphere boundary
class OuterBoundary(SubDomain):
    def inside(self, x, on_boundary):
        r2 = sqrt(x[0]**2 + x[1]**2 + x[2]**2)
        print r2
        return abs(r2-oradius) < tol

#Define the inner sphere boundary
class InnerBoundary(SubDomain):
    def inside(self, x, on_boundary):
        r2 = sqrt(x[0]**2 + x[1]**2 + x[2]**2)
        print r2
        return abs(r2-iradius) < tol

boundary_parts = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundary_parts.set_all(0)
outer_boundary = OuterBoundary()
inner_boundary = InnerBoundary()

outer_boundary.mark(boundary_parts, 1)
outer_boundary.mark(boundary_parts, 2)
ds = Measure("ds")[boundary_parts]

print np.linalg.norm(boundary_parts.array())
plot(boundary_parts)
interactive()
#Function Space
VV = VectorFunctionSpace(mesh, "Lagrange", 2, 3) # displacement space

u = TrialFunction(VV)
v = TestFunction(VV)

eps_u = sym(grad(u))
eps_v = sym(grad(v))
sig_u = lmbda*tr(eps_u)*Identity(3) + 2.*mu*eps_u

# Traction vector [1., 1., 1.] applied on the outer boundary
t = Constant((1., 1., 1.))

F = inner(sig_u, eps_v)*dx - inner(t, v)*ds(1) 
a, L = lhs(F), rhs(F)

u_h = Function(VV)

A, b = assemble_system(a, L)

print b.norm('l2')

solve(A, u_h.vector(), b, 'lu')

ux, uy, uz = u_h.split()

p
