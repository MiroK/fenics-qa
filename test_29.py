from dolfin import *
import sys, math, numpy

mesh = UnitIntervalMesh(100)

subdomains = FacetFunction('size_t', mesh, 1) #MeshFunction should be used?

class domain_1(SubDomain):
    def inside(self, x, on_boundary):
        return True if 0.3 <= x[0] <= 0.5 else False

subdomains.set_all(0)
subdomain1 = domain_1()
subdomain1.mark(subdomains, 1)

dx=Measure('dx')[subdomains]
ds=Measure('ds')[subdomains]

V=FunctionSpace(mesh,"CG",1)

u=TrialFunction(V)
v=TestFunction(V)

f1=Expression('x[0]',cell=interval)
f2=Expression('3*x[0]',cell=interval)

s=Function(V)

n=FacetNormal(mesh)

a=dot(grad(u),grad(v))*dx(0)+dot(grad(u),grad(v))*dx(1)
L=dot(inner(grad(f1),n),v)*ds(0)+dot(inner(grad(f2),n),v)*ds(1)

tol = 1E-16
def left_boundary(x, on_boundary):
    return on_boundary and abs(x[0]) < tol

def right_boundary(x, on_boundary):
    return on_boundary and abs(x[0]-1) < tol

Gamma_0_delr0 = DirichletBC(V, Constant(0.0), left_boundary)
Gamma_1_delr0 = DirichletBC(V, Constant(2.6), right_boundary)
bcs_delr0 = [Gamma_0_delr0, Gamma_1_delr0] #don't know how to impose the

solve(a==L,s,bcs_delr0)

plot(s)
interactive()
