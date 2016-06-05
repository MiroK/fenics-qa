from dolfin import *

import numpy, sys

nx = 10
degree = 1
mesh = UnitSquareMesh(nx, nx)

V2 = FunctionSpace(mesh, "DG",degree)
V0 = FunctionSpace(mesh, "CG", degree)
W = MixedFunctionSpace([V0, V2,V2])

tol = DOLFIN_EPS
def left_boundary(x, on_boundary):
    return on_boundary and abs(x[0]) < tol
def low_boundary(x, on_boundary):
    return on_boundary and abs(x[1]) < tol
def right_boundary(x, on_boundary):
    return on_boundary and abs(x[0]-1) < tol
def up_boundary(x, on_boundary):
    return on_boundary and abs(x[1]-1) < tol

pex = Expression('x[0]*x[1]*(1.0 -x[0])*(1.0 -x[1]) ')
Gamma_0 = DirichletBC(W.sub(0), pex , left_boundary)
Gamma_1 = DirichletBC(W.sub(0), pex , right_boundary)
Gamma_2 = DirichletBC(W.sub(0), pex, up_boundary)
Gamma_3 = DirichletBC(W.sub(0), pex , low_boundary)

bcs = [Gamma_0, Gamma_1, Gamma_2, Gamma_3]

(u, q1, q2) = TrialFunctions(W)
(ut, q1t, q2t) = TestFunctions(W)
w = Function(W)
u, q1,q2 = ( w[0], w[1], w[2])

alfa = 1.E-10
beta = 1.E-10

h = Expression('(2.0 -x[0]) ')
uini = Expression(' x[0]*x[1]*(1.0 -x[0])*(1.0 -x[1]) ')

px = Expression('x[1]*(1.0 -x[1])*(1.0 -2.0*x[0]) ')
fsec = -pex + h*px

wini = Function(W)
q1ex = Expression('x[1]*(1.0 -x[1])*(1.0 -2.0*x[0]) ')
q2ex = Expression('x[0]*(1.0 -x[0])*(1.0 -2.0*x[1]) ')
wini0 = interpolate(uini,V0)
wini1 = interpolate(q1ex,V0)
wini2 = interpolate(q2ex,V0)

a1 = -inner(h*u , Dx(ut,0))*dx - inner( u,Dx(ut,1))*dx
a2 = inner(q1,q1t)*dx + inner(-alfa*h*u,Dx(q1t,0))*dx
a3 = inner(-alfa*h*u,Dx(q2t,1))*dx + inner(q2,q2t)*dx
L = fsec*ut*dx
F = a1 - L + a2 + a3

J = derivative(F,w)
wini = interpolate(Expression(('w0', 'w1', 'w2'), w0=wini0, w1=wini1, w2=wini2),
                   W)

problem = NonlinearVariationalProblem(F, wini, bcs=bcs,J=J)
solver = NonlinearVariationalSolver(problem)
solver.solve()
