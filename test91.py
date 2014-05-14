from dolfin import *

mesh = CircleMesh(Point(0., 0.), 1, 0.1)

# params
dt = Constant(0.1)  # time step
mu = Constant(100)   # diffusion
sigma = 0.5          # magnitude of delta function

# refine mesh aroung the support of sourse
for ref_level in range(3):
  cell_f = CellFunction('bool', mesh, False)
  for cell in cells(mesh):
    x = cell.midpoint().x()
    y = cell.midpoint().y()
    r = x**2 + y**2
    if r < 0.1:
      cell_f[cell] = True
  mesh = refine(mesh, cell_f)

# approximation to delta function
class Delta(Expression):
  def __init__(self, a):
    self.a = a

  def eval(self, values, x):
    a = self.a
    r = x[0]**2 + x[1]**2
    values[0] = exp(-r/a**2)/sqrt(pi)/a

# var. formulation
V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

u0 = interpolate(Constant(0.), V)  # initial conditions
bc = DirichletBC(V, Constant(0.), DomainBoundary())

f = Delta(1./10)
plot(f, interactive=True, mesh=mesh, title='source')

F = (inner(u - u0, v)/dt)*dx + mu*inner(grad(u), grad(v))*dx - inner(f, v)*dx
a = lhs(F)
L = rhs(F)

u = Function(V)
t = 0
while t < 1:
  t += float(dt)
  A, b = assemble_system(a, L, bc)

  solve(A, u.vector(), b)

  u0.assign(u)
  plot(u0, title='solution', interactive=True)

  print t, u0.vector().max()
