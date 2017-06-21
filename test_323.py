from dolfin import *

mesh = UnitSquareMesh(6,4)

ele = FiniteElement("CG", mesh.ufl_cell(), 1)

V = FunctionSpace(mesh, MixedElement([ele,ele]))

u, v = TrialFunction(V), TestFunction(V)

u_r, u_i = split(u)

v_r, v_i = split(v)

f_r = Constant(0.)

f_i = Constant(0.)

def u0_boundary(x, on_boundary):

    return on_boundary

bc = DirichletBC(V, Constant((0,0)), u0_boundary)

a_r = inner(grad(u_r),grad(v_r)) * dx - inner(grad(u_i),grad(v_i)) * dx

a_i = inner(grad(u_r),grad(v_i)) * dx + inner(grad(u_i),grad(v_r)) * dx

L_r = inner(f_r,v_r) * dx - inner(f_i,v_i) * dx

L_i = inner(f_r, v_i) * dx + inner(f_i,v_r) * dx

a = a_r+a_i

L = L_r+L_i

A = assemble(a)

b = assemble(L)

bc.apply(A,b)

ps_r = PointSource(V, Point(.5,.5))
#ps_i = PointSource(V.sub(1), Point(.5,.5))

# ps_i.apply(b)
ps_r.apply(b)


print b.array() # inspect ps has been applied

u = Function(V)

solve(A, u.vector(), b, 'lu')

u_r, u_i = u.split(True)

plot(u_i, title='imag');plot(u_r);interactive()
