from dolfin import *

mesh = UnitSquareMesh(4, 4)

T = TensorElement('Lagrange', mesh.ufl_cell(), 1)
T = FunctionSpace(mesh, T)

sigma = interpolate(Constant(((1, 2), (4, 5))), T)

e1 = Constant((1, 0))
e2 = Constant((0, 1))

s12 = e1[i]*sigma[i, j]*e2[j] #inner(e1, dot(sigma, e2))
s21 = e2[i]*sigma[i, j]*e1[j] #inner(e2, dot(sigma, e1))
f = as_vector((s12, s21))

V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

a = -div(f*u)*v*dx
A = assemble(a)

print '|A|', A.norm('linf')
