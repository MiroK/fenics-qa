from dolfin import *

# ------------------------
# BEGIN: class definitions

class LowerBoundaryOfUnitSquare(SubDomain):
    """ This class represents and manipulates 
    the lower boundary of the unit square
    [0,1]x[0,1]. """
    def inside(self, x, on_boundary):
        tol = 1E-14   # tolerance for coordinate comparisons
        return on_boundary and abs(x[1]) < tol

class UpperBoundaryOfUnitSquare(SubDomain):
    """ This class represents and manipulates 
        the upper boundary of the unit square
        [0,1]x[0,1]. """
    def inside(self, x, on_boundary):
        tol = 1E-14   # tolerance for coordinate comparisons
        return on_boundary and abs(x[1] - 1) < tol

class LeftBoundaryOfUnitSquare(SubDomain):
    """ This class represents and manipulates 
        the left boundary of the unit square
        [0,1]x[0,1]. """
    def inside(self, x, on_boundary):
        tol = 1E-14   # tolerance for coordinate comparisons
        return on_boundary and abs(x[0]) < tol

class RightBoundaryOfUnitSquare(SubDomain):
    """ This class represents and manipulates 
        the right boundary of the unit square
        [0,1]x[0,1]. """
    def inside(self, x, on_boundary):
        tol = 1E-14   # tolerance for coordinate comparisons
        return on_boundary and abs(x[0] - 1) < tol

# END: class definitions
# ----------------------

# --------------------
# BEGIN: Main function

nx = 64
ny = 64
mesh = UnitSquareMesh(nx, ny)

# create a mesh function over cell facets
mesh_fn_marking_bdry_parts = MeshFunction("size_t", mesh, mesh.topology().dim()-1)

Gamma_0 = LowerBoundaryOfUnitSquare()
Gamma_0.mark(mesh_fn_marking_bdry_parts, 0)
Gamma_1 = RightBoundaryOfUnitSquare()
Gamma_1.mark(mesh_fn_marking_bdry_parts, 1)
Gamma_2 = UpperBoundaryOfUnitSquare()
Gamma_2.mark(mesh_fn_marking_bdry_parts, 2)
Gamma_3 = LeftBoundaryOfUnitSquare()
Gamma_3.mark(mesh_fn_marking_bdry_parts, 3)

# set up and solve the weak formulation of the Neumann 
# boundary-value problem for the Laplacian
V = FunctionSpace(mesh, "CG", 1)
R = FunctionSpace(mesh, "R", 0)
W = V * R

(u, c) = TrialFunction(W)
(v, d) = TestFunction(W)

f = Constant("0")
g_0 = Constant("0")
g_1 = Constant("1")
g_2 = Constant("0")
g_3 = Constant("-1")

(u, c) = TrialFunction(W)
(v, d) = TestFunction(W)

a = ( inner(grad(u),grad(v)) + c*v + d*u )*dx
L = f*v*dx + g_0*v*ds(0, subdomain_data=mesh_fn_marking_bdry_parts) \
           + g_1*v*ds(1, subdomain_data=mesh_fn_marking_bdry_parts) \
           + g_2*v*ds(2, subdomain_data=mesh_fn_marking_bdry_parts) \
           + g_3*v*ds(3, subdomain_data=mesh_fn_marking_bdry_parts) 

A = assemble(a)
b = assemble(L)

w = Function(W)
solve(A, w.vector(), b, 'lu')
(u, c) = w.split()

# plot the numerically obtained solution
plot(u, title = 'The solution u')
plot(Expression('x[0]-0.5'), mesh=mesh)

plot(grad(u))

interactive()

# compute and print the flux through the boundary
unitNormal = FacetNormal(mesh)

for i in range(0,4):
    flux = dot(grad(u), unitNormal)*ds(i,subdomain_data=mesh_fn_marking_bdry_parts)

    flux = assemble(flux)
    print "outward flux of u through Gamma_{0}: ".format(i), flux

# END: Main function
# --------------------------
