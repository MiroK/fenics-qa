from dolfin import *
mesh = UnitSquareMesh(32,32)
#Function Spaces
V = VectorFunctionSpace(mesh, "Lagrange", 2)
Q = FunctionSpace(mesh, "Lagrange", 1)
W = V*Q

#Trial and Test Functions
(u,p) = TrialFunctions(W)
(v,q) = TestFunctions(W)

#Exact Values / RHS values

u_ex = Expression(('cos(pi*x[1])', 'sin(pi*x[0])'))
p_ex = Expression('pi*cos(pi*x[0])*cos(pi*x[1])')

f = Expression(('-pi*pi*sin(pi*x[0])*cos(pi*x[1]) + pi*pi*cos(pi*x[1])',
                'pi*pi*sin(pi*x[0]) - pi*pi*sin(pi*x[1])*cos(pi*x[0])'))

#Boundary Conditions
def right(x, on_boundary): return x[0] > 1.0 - DOLFIN_EPS
def left(x, on_boundary): return x[0] < DOLFIN_EPS
def top(x, on_boundary): return x[1] > 1.0 - DOLFIN_EPS
def bottom(x, on_boundary): return x[1] < DOLFIN_EPS

bc_u_right = DirichletBC(W.sub(0), u_ex, right)
bc_u_left = DirichletBC(W.sub(0), u_ex, left)
bc_u_top = DirichletBC(W.sub(0), u_ex, top)
bc_u_bot = DirichletBC(W.sub(0), u_ex, bottom)

bcs = [bc_u_right, bc_u_left, bc_u_top, bc_u_bot]

#Variational Problem
a = inner(nabla_grad(u), nabla_grad(v))*dx - div(v)*p*dx + q*div(u)*dx
L = inner(f, v)*dx

#Solve
w = Function(W)
solve(a == L, w, bcs)
(u, p) = w.split(True)
normalize(p.vector())

#Plot
plot(u, title="My Velocity")
plot(u_ex,mesh, title="Exact Velocity")

plot(p, title="My Pressure")
plot(p_ex, mesh, title="Exact Pressure")

interactive()
