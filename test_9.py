import dolfin
from mshr import *
def set_mesh(ref,h,w):
   p1 = dolfin.Point(0.,0.)
   p2 = dolfin.Point(w,h)
   p3 = dolfin.Point(w,1/2.0*h)
   domain =   Rectangle(p1,p2)
   domain.set_subdomain(1, Rectangle(p1, p3))
   mesh2d = generate_mesh(domain, ref)
   return (mesh2d)
   

from dolfin import *

mesh= set_mesh(10,1.,1.2)

V = VectorFunctionSpace(mesh, "Lagrange", 1)

Hi, it is not clear what your problem is but I suspect you get warnings about `bc2` being unable to find facets to
apply boundary conditions at. This is caused by the definitions of your subdomains where you use `==` to compare two
doubles. You should use instead the `near` function which checks equality of its arguments within some tolerance. That
is your domains could look as follows 

TOL = DOLFIN_EPS

class Corner(SubDomain):
  def inside(self, x, on_boundary):
    return near(x[0], 0., TOL) and near(x[1], 0., TOL)

class Border(SubDomain):
  def inside(self, x, on_boundary):
    return near(x[1], 0., TOL)

border = Border()
corner = Corner()
bc1 = DirichletBC(V, Constant((0.,0.)),corner,method="pointwise")
bc2 = DirichletBC(V.sub(1), 0., border)
bcs = [bc1, bc2]

u = TrialFunction(V)
v = TestFunction(V)

M = assemble(inner(u, v)*dx)
bc1.apply(M)
bc2.apply(M)
