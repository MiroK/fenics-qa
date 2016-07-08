from dolfin import *

mesh = UnitCubeMesh(2, 2, 2)
V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)
lamb = 1.

x = SpatialCoordinate(mesh)
u_exact = (sin(x[2]))**lamb*sin(lamb*x[1]) 
f = lamb*(lamb + 1.)*(sin(x[2]))**lamb*sin(lamb*x[1])
a = ( Dx(u,0)*Dx(v,0) + 1./(x[0])**2*Dx(u,1)*Dx(v,1) \
             + 1./(x[0]*sin(x[1]))**2*Dx(u,2)*Dx(v,2) )*dx
l = f*v*dx
A = assemble(a)
L = assemble(l)
