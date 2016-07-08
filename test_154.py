from dolfin import *

mesh = UnitCubeMesh(3, 3, 3)
V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

r = Expression("x[0]", degree=1)
phi = Expression("x[1]", degree=1)
theta = Expression("x[2]", degree=1)
lamb = 2
f = Expression("lamb*(lamb + 1)*(sin(x[2]))*lamb*sin(lamb*x[1])", lamb=lamb)  
a = ( Dx(u,0)*Dx(v,0) + (Constant(1.)/r)**2*Dx(u,1)*Dx(v,1) +     Constant(1.)/(r*sin(phi))**2*Dx(u,2)*Dx(v,2) )*dx
L = inner(f,v)*dx
