from mshr import Sphere, generate_mesh
from dolfin import *

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

#Function Space
VV = VectorFunctionSpace(mesh, "Lagrange", 1, 3) # displacement space
Q = VectorFunctionSpace(mesh, 'R', 0, 6)
W = MixedFunctionSpace([VV, Q])

u, p = TrialFunctions(W)
v, q = TestFunctions(W)

eps_u = sym(grad(u))
eps_v = sym(grad(v))
sig_u = lmbda*tr(eps_u)*Identity(3) + 2.*mu*eps_u

t = Constant((1, 1, 1))

# Nullspace of rigid motions
# Translation
Z_transl = [Constant((1, 0, 0)), Constant((0, 1, 0)), Constant((0, 0, 1))]
# Rotations
Z_rot = [Expression(('0', 'x[2]', '-x[1]')),
         Expression(('-x[2]', '0', 'x[0]')),
         Expression(('x[1]', '-x[0]', '0'))]
# All
Z = Z_transl + Z_rot

a = inner(sig_u, eps_v)*dx\
   -sum(p[i]*inner(v, Z[i])*dx for i in range(len(Z)))\
   -sum(q[i]*inner(u, Z[i])*dx for i in range(len(Z))) 

L = inner(t, v)*ds

w_h = Function(W)

A, b = assemble_system(a, L)
print b.norm('l2')

solve(A, w_h.vector(), b, 'lu')

uh, ph = w_h.split(deepcopy=True)

plot(uh, mode = "displacement", interactive=True, wire_frame=True, axes=True,rescale=False,scalarbar=True)
