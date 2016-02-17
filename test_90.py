from dolfin import *

mesh = UnitCubeMesh(3, 3, 3)
V = VectorFunctionSpace(mesh,'CG', 1)
dA = Measure('ds')
dV = Measure('dx')
p0 = 101325.0E-6# in MPa
g = 9.81# in m/s2
rho_oel = 860# in kg/(m)3
h = 1.5350# in m
padd = Expression(('A*B*(C-x[2]*1.0E-3)'),A=rho_oel, B=g, C=h) # in 1.0E-6MPa
pin = p0 + 1.0E-6*padd# in MPa
n = FacetNormal(mesh)
tr_a = -p0*n
tr_i = -pin*n
null = Constant((1.0, 0.0, -10))
bc = [DirichletBC(V, null, 'on_boundary')]
u = TrialFunction(V)
del_u = TestFunction(V)
nu = 0.38
E = 1400.0 # MPa
G = E/(2.0*(1.0+nu))
l = 2.0*G*nu / (1.0 - 2.0*nu)
mu = G
delta = Identity(3)

eps = lambda u: as_tensor(1.0/2.0*(u[i].dx(j) + u[j].dx(i)), (i,j))
sigma = lambda u: as_tensor(l*eps(u)[k,k]*delta[i,j] + 2.0*mu*eps(u)[i,j], (i,j))
a = sigma(u)[j,i]*del_u[i].dx(j)*dV
L = tr_a[i]*del_u[i]*dA + tr_i[i]*del_u[i]*dA

disp = Function(V)
solve(a == L, disp, bcs=bc)


sdev = as_tensor(sigma(disp)[i,j]- 1.0/3.0*sigma(disp)[k,k]*delta[i,j], (i,j))
mises = as_tensor((3.0/2.0*sdev[i,j] * sdev[i,j])**0.5, () )
seq = project(mises, FunctionSpace(mesh, 'CG', 1), 
              solver_type='mumps', 
              form_compiler_parameters={"cpp_optimize":True,"representation":"quadrature","quadrature_degree":2})
print seq.vector().norm('l2')
