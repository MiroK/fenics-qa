from dolfin import *
from scipy.linalg import eigvals

mesh = UnitTriangleMesh()
#mesh = UnitSquareMesh(1, 1)

Vu = VectorFunctionSpace(mesh, 'DG', 1)
Vv = FunctionSpace(mesh, 'DGT', 1) 
u = TrialFunction(Vu)
v = TestFunction(Vv)

n = FacetNormal(mesh)
w = Constant([1, 1])

a = dot(u, n)*v*ds
L = dot(w, n)*v*ds
for R in '+-':
    a += dot(u(R), n(R))*v(R)*dS
    L += dot(w(R), n(R))*v(R)*dS

A = assemble(a)
b = assemble(L)
U = Function(Vu)

if True:
    ls = LocalSolver(a, L)
    ls.factorize()

    print A.array()
    print eigvals(A.array())

    ls.solve_local_rhs(U)
else:
    solve(A, U.vector(), b)
    print U.vector().array()
