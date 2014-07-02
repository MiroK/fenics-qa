from dolfin import *
from numpy.linalg import norm as np_norm

# Mesh and functionspace
mesh = UnitSquareMesh(100,100)
V = VectorFunctionSpace(mesh, "CG", 1)
fe = Function(V)
ue = Function(V)

# Extract x and y coordinates of mesh and
# align with dof structure
dim = V.dim()
N = mesh.geometry().dim()
coor = V.dofmap().tabulate_all_coordinates(mesh).reshape(dim,N)
fx_dofs = V.sub(0).dofmap().dofs()
fy_dofs = V.sub(1).dofmap().dofs()
x = coor[:, 0]   # x for fx and fy
y = coor[:, 1]   # y for fx and fy
fx_x, fx_y = x[fx_dofs], y[fx_dofs]  # x, y of components
fy_x, fy_y = x[fy_dofs], y[fy_dofs]

print V.dim()/2, len(fx_x)
# x and y components of vector function
fx = fx_x*fx_y
fy = 100+0*fy_x

# Insert values of fx and fy into the function fe
fe.vector()[fx_dofs] = fx
fe.vector()[fy_dofs] = fy


# Function in fenics code for testing purposes
func = Expression(("x[0]*x[1]","100"))
f = interpolate(func, V)
ue.assign(f)

# Check that the methods give the same result
ufunc = fe.vector().array()
uexpr = ue.vector().array()

print ufunc
print uexpr
print ufunc - uexpr
print 'Match?', near(np_norm(ufunc-uexpr), DOLFIN_EPS)
