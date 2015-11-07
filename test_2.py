from dolfin import *
mesh = IntervalMesh(20, 0, 1)

def DirichletBoundary(x, on_boundary):
    return on_boundary

u0=Expression('0.1*x[0]')
v0=Expression('0.0')

FS1 = FunctionSpace(mesh, "CG", 4)
FS2 = FunctionSpace(mesh, "CG", 4)
W = FS1 * FS2

bc = [DirichletBC(W.sub(0), u0, DirichletBoundary),DirichletBC(W.sub(1), v0, DirichletBoundary)]

(u,v) = TrialFunctions(W)
(du,dv) = TestFunctions(W)

K1=Expression('100.0*(x[0]+1.0)')
K2=Expression('500.0*(x[0]+1.0)')
dK1=Expression('100.0')
dK2=Expression('500.0')
ddK1=Expression('0.0')
ddK2=Expression('0.0')

a= (
    -inner(grad(u),grad(du))
    -inner(v,du)
    #~ -2.0*dK2*v*grad(v)[0]
    -K2*inner(grad(v),grad(dv))
    +ddK2*inner(v,dv)
    #~ -dK1*grad(u)[0]*dv
    -K1*inner(v,dv)
    )*dx

f = Constant(0.0)
L = f*dv*dx


w = Function(W)
solve( a == L, w, bc)

(u, v) = w.split()

plot(u, interactive = True)
