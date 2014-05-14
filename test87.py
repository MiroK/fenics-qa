from dolfin import *
from math import sin, cos

mesh = RectangleMesh(-1, -1, 1, 1, 100, 100)
V = FunctionSpace(mesh, "CG", 1)
u = TrialFunction(V)
v = TestFunction(V)

a = inner(grad(u), grad(v))*dx
L = Constant(0)*v*dx
bc = DirichletBC(V, Constant(0), DomainBoundary())
A, b = assemble_system(a, L, bc)

u = Function(V)
n_steps = 20
dt = 1./n_steps
t = 0
for step in range(n_steps):
    cn = cos(2*pi*step*dt)
    sn = sin(2*pi*step*dt)
    den = 1 + sn**2

    x = 0.5*cn/den
    y = 0.5*cn*sn/den
    t += dt

    delta = PointSource(V, Point(x, y), abs(x*y)*t*100)
    delta.apply(b)
    solve(A, u.vector(), b)

    plot(u, interactive=True)

    # reset rhs
    b.zero()
    bc.apply(b)
