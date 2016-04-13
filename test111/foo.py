from  dolfin  import *

mesh = Mesh('moebius.xml.gz')
V = FunctionSpace(mesh, "Lagrange", 1)

print V.element().geometric_dimension(), V.element().topological_dimension()
print V.ufl_cell()

u0 = Constant(0.0)
bc = DirichletBC(V, u0, "on_boundary")

u = TrialFunction(V)
v = TestFunction(V)
f = Constant(1.0)
a = inner(grad(u), grad(v))*dx
L = f*v*dx

u = Function(V)
solve(a == L, u, bc)

plot(u, interactive=True)
