from dolfin import *

mesh = UnitCubeMesh(3, 3, 3)
V  = VectorFunctionSpace(mesh, "Lagrange", 1)  

r1 = 4.0/5.0-DOLFIN_EPS;

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0)

left = Left()

if x[0] >= r1: 
    f= Expression((0.0, 0.0, -10*exp(2*x[2]+x[1])))
else:
    f= Expression((0.0, 0.0, 0.0))

E = 200*(10**9)
nu = 0.3

mu, lmbda = E/(2.0*(1.0 + nu)), E*nu/((1.0 + nu)*(1.0 - 2.0*nu))

sigma = 2*mu*sym(grad(u)) + lmbda*tr(grad(u))*Identity(v.cell().d)   

a = inner(sigma, grad(v))*dx 
L = dot(f, v)*dx                          # body force f   

c = Constant((0.0, 0.0, 0.0))
bc = [DirichletBC(V, c, left)]       #gmsh tag for fixed face is 1 but doing it the conventional way

