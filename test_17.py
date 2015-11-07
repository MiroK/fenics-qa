from dolfin import *
mesh = UnitSquareMesh(8,8)

V = FunctionSpace(mesh,"Lagrange",1)
ME = MixedFunctionSpace([V,V,V])
u = Function(ME)
c, mu, n = u.split() 

z = Function(V)
z.interpolate(Expression("sin(x[0]*x[1])"))

# to ME0 from V
assigner0 = FunctionAssigner(ME.sub(0), V)
# assign to c function z
assigner0.assign(c, z)

plot(z, title='z')
plot(c, title='c', interactive=True)

u.vector()[:] *= 10

assigner1 = FunctionAssigner(V, ME.sub(0))
assigner1.assign(z, c)

plot(z, title='z')
plot(c, title='c', interactive=True)
