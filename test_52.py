from dolfin import *
import sys, math, numpy

nx=10
ny=20
nz=20

mesh = UnitCubeMesh(nx,ny,nz)
V= FunctionSpace(mesh, 'Lagrange', 1)

class LeftBoundary(SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and near(x[0],0.0)

class RightBoundary(SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and near(x[0],1.0)

left_bndr = LeftBoundary()
right_bndr = RightBoundary()

Gamma_L = DirichletBC(V, Constant(0),left_bndr)
Gamma_R = DirichletBC(V, Constant(1.0), right_bndr)

u = TrialFunction(V)
v = TestFunction(V)
f= Constant(0)

a = inner(nabla_grad(u),nabla_grad(v))*dx
L= f*v*dx
u= Function(V)
solve(a == L, u, [Gamma_L, Gamma_R])

boundary_parts = FacetFunction("size_t", mesh)
left_bndr.mark(boundary_parts,1)
right_bndr.mark(boundary_parts,2)

class InnerBoundary(SubDomain):
    def inside(self,x, on_boundary):
        return near(x[0],0.5)

inner_bndr = InnerBoundary()
inner_bndr.mark(boundary_parts,3)
ds = Measure("ds")[boundary_parts]
dS = Measure("dS")[boundary_parts]

flux_l=dot(Constant((1,0,0)),nabla_grad(u))*ds(1)
flux_r=dot(Constant((1,0,0)),nabla_grad(u))*ds(2)
print "flux left: ", assemble(flux_l)
print "flux right: ", assemble(flux_r)

flux_c = dot(Constant((1.0,0.0,0)), avg(nabla_grad(u)))*dS(3)
print "Central flux: ", assemble(flux_c)
