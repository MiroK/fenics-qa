from dolfin import *

mesh = RectangleMesh(Point(0,0),Point(1, 1),40,20)
E0 = 10e5
nu = 0.28
alpha = 5e-6
Tref = 293
T = 1273
delta_T = T - Tref
E  = E0
mu = E/(2*(1+nu))
lamb = E*nu/((1+nu)*(1-2*nu))
thermal_coef= alpha*delta_T*(3+2*mu/lamb)
# define rectangular beam
# define function space
V = VectorFunctionSpace(mesh,'CG',1)
# Boundary Conditions
tol = 1e-14
def clamped_boundary(x,on_boundary):
    return on_boundary and x[0]<tol
bc = DirichletBC(V,Constant((0,0)),clamped_boundary)
# define strain-displacement 
def epsilon(u):
    return 0.5*(nabla_grad(u)+nabla_grad(u).T)
# define stifness tensor
def sigma(u):
    return lamb*((nabla_div(u)-alpha*delta_T*(3+2*mu/lamb))*Identity(d))+2*mu*epsilon(u)
# elasticity problem
u = TrialFunction(V)
d = u.geometric_dimension()
v = TestFunction(V)
f = Expression(( '0','-20'))
form = inner(sigma(u),epsilon(v))*dx - dot(f, v)*dx
a = lhs(form)
L = rhs(form)

u = Function(V)
solve(a == L,u,bc)
plot(u,mode='displacement')
