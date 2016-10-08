from dolfin import *
mesh = Mesh('mesh_spherical.xml')    

V = FunctionSpace(mesh, "CG", 1)
u = TrialFunction(V)
v = TestFunction(V)

lamb = .6
x = SpatialCoordinate(mesh)
f = lamb*(lamb + 1.)*(sin(x[2]))**lamb*sin(lamb*x[1])
a = Dx(u,0)*Dx(v,0)*dx
a += 1./(x[0])**2*Dx(u,1)*Dx(v,1)*dx
a += 1./(x[0]*sin(x[1]))**2*Dx(u,2)*Dx(v,2)*dx
l = f*v*dx
A = assemble(a)
#L = assemble(l)
