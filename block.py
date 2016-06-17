from dolfin import *

# Load mesh and subdomains
mesh = Mesh("dolfin_fine.xml.gz")
sub_domains = MeshFunction("size_t", mesh, "dolfin_fine_subdomains.xml.gz")

# Define function spaces
V1 = VectorFunctionSpace(mesh, 'CG', 1)
V2 = VectorFunctionSpace(mesh, 'Bubble' ,3)
Q = FunctionSpace(mesh, 'CG', 1)

u, b, p = map(TrialFunction, (V1, V2, Q))
v, c, q = map(TestFunction, (V1, V2, Q))

# No-slip boundary condition for velocity
# x1 = 0, x1 = 1 and around the dolphin
noslip = Constant((0, 0))
bc0 = DirichletBC(V1, noslip, sub_domains, 0)

# Inflow boundary condition for velocity
# x0 = 1
inflow = Expression(("-sin(x[1]*pi)", "0.0"), degree=2)
bc1 = DirichletBC(V1, inflow, sub_domains, 1)

# Upper triangular part of A
a00 = inner(grad(u), grad(v))*dx
a01 = inner(grad(b), grad(v))*dx
a02 = inner(p, div(v))*dx
a11 = inner(grad(b), grad(c))*dx
a12 = inner(p, div(c))*dx

f = Constant((0, 0))
# Rhs components
L0 = inner(f, v)*dx
L1 = inner(f, c)*dx

from block import block_assemble

bb = block_assemble([L0, L1, 0])


# # Collect boundary conditions
# bcs = [bc0, bc1]
# 
# # Define variational problem
# (u, p) = TrialFunctions(W)
# (v, q) = TestFunctions(W)
# f = Constant((0, 0))
# a = (inner(grad(u), grad(v)) - div(v)*p + q*div(u))*dx
# L = inner(f, v)*dx
# 
# # Compute solution
# w = Function(W)
# solve(a == L, w, bcs)
# 
# # Split the mixed solution using deepcopy
# # (needed for further computation on coefficient vector)
# (u, p) = w.split(True)
# 
# print("Norm of velocity coefficient vector: %.15g" % u.vector().norm("l2"))
# print("Norm of pressure coefficient vector: %.15g" % p.vector().norm("l2"))
# 
# # # Split the mixed solution using a shallow copy
# (u, p) = w.split()
# 
# # Save solution in VTK format
# ufile_pvd = File("velocity.pvd")
# ufile_pvd << u
# pfile_pvd = File("pressure.pvd")
# pfile_pvd << p
# 
# # Plot solution
# plot(u)
# plot(p)
# interactive()
