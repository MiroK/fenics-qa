from dolfin import *
import numpy as np

# Create mesh and define function space
mesh = UnitSquareMesh(4, 4)
V = FunctionSpace(mesh, "Lagrange", 1)

u_e=Expression('1+x[0]*x[0]+2*x[1]*x[1]')	                 #exact solutin	

# Define Dirichlet boundary (x = 0 or x = 1)
class Left(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0]*(1-x[0]),0.0)

#Define the right dirichlet boundary condition
class Right(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1]*(1-x[1]), 0.0)

left=Left()
right=Right()


# Define boundary condition
u0 = Expression('1+x[0]*x[0]+2*x[1]*x[1]')	
bc = DirichletBC(V, u0, left)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression("-6")	
ur = Expression('4*x[1]')	
a = inner(grad(u), grad(v))*dx
L = f*v*dx + ur*v*ds

# Compute solution
u = Function(V)                            # u is the solution with CG method  
solve(a == L, u, bc)


u_e_Ve = interpolate(u_e, V)

error = (u - u_e_Ve)**2*dx
k=sqrt(assemble(u_e_Ve**2*dx))
E = assemble(error)
print E
k=sqrt(assemble(u_e_Ve**2*dx))         #to get relative L2 norm

#print k
#print E
print('L2 norm using CG Method : ',E/k)

#plot(u)
#plot(u_e_Ve)
#interactive()
