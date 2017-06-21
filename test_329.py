from dolfin import *
import mshr
domain = mshr.Rectangle(Point(0,0), Point(1,1))-mshr.Circle(Point(0.5,0.5), 0.1)
mesh = mshr.generate_mesh(domain, 10)
V = VectorFunctionSpace(mesh, 'CG',1)
u = TrialFunction(V)
v = TestFunction(V)
lmbda=121.15E3
mu=80.77E3
def epsilon(u):
    return 0.5*(grad(u)  + grad(u).T)

def sigma(u):
    return 2.0*mu*epsilon(u) + lmbda*tr(epsilon(u))*Identity(2)

class top(SubDomain):
    def inside(self,x,on_boundary):
        tol = 1e-5
        return abs(x[1]-1.0) < tol and on_boundary

class bottom(SubDomain):
    def inside(self,x,on_boundary):
        tol = 1e-5
        return abs(x[1]) < tol and on_boundary
Bottom = bottom()
Top = top()
u0=Constant(("0.0","0.0"))
bc1=DirichletBC(V,u0,Bottom)
bc2=DirichletBC(V.sub(1),Constant(1.0),Top)
bc = [bc1, bc2]

Eu = inner(nabla_grad(v), sigma(u))*dx
uh = Function(V)

M = uh[0]*dx
tol = 1.e-3

solver_parameters = {"error_control":{"dual_variational_solver":{"linear_solver": "cg"}}}
solve(lhs(Eu) == rhs(Eu), uh, bc , tol=tol , M=M,solver_parameters=solver_parameters)

plot(u.root_node (), title="Solution  on  initial  mesh")
plot(u.leaf_node (), title="Solution  on final  mesh")
interactive ()
