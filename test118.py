from dolfin import *

class Dielectric(SubDomain):
    def inside(self, x, on_boundary):
        return (between(x[1], (0.0, 0.25)))
dielectric = Dielectric()

mesh = RectangleMesh ( 0, 0, 1, 0.5, 8, 4)
domains = CellFunction("size_t", mesh)
domains.set_all(0)
dielectric.mark(domains, 1)
dx = Measure("dx")[domains]

V = FunctionSpace ( mesh, "Nedelec 1st kind H(curl)", 3 )
u = TestFunction(V)
v = TrialFunction(V)

def curl_t(w):
    return Dx(w[1], 0) - Dx(w[0], 1)

e1 = Constant(0)
e2 = Constant(4)

s = curl_t(u)*curl_t(v)*dx(0) + curl_t(u)*curl_t(v)*dx(1)
t = dot( v, u ) * ( 1 - e1) * dx(0) + dot( v, u) * ( 1 - e2) * dx(1)

S = PETScMatrix()
T = PETScMatrix()

class ElectricWall(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary;

assemble(s, tensor=S)
assemble(t, tensor=T)

electric_wall = DirichletBC( V, Expression(("0.0", "0.0")), ElectricWall())
electric_wall.apply (S)
electric_wall.apply (T)
























