from dolfin import *
from slepc4py import SLEPc

mesh = UnitSquareMesh(40, 40)
V = FunctionSpace(mesh, 'Lagrange', 2)
u = TrialFunction(V)
v = TestFunction(V)
bc = DirichletBC(V, Constant(0.), 'on_boundary')

a = inner(grad(u), grad(v))*dx
m = inner(u, v)*dx
L = inner(Constant(0), v)*dx

A, M = PETScMatrix(), PETScMatrix()
b = PETScVector()
assemble_system(a, L, A_tensor=A, b_tensor=b)
assemble_system(m, L, A_tensor=A, b_tensor=b)

E = SLEPc.EPS().create()
E.setOperators(A.mat(), M.mat())
E.setInterval(0.5, 2)
E.setProblemType(SLEPc.EPS.ProblemType.HEP)
E.setWhichEigenpairs(SLEPc.EPS.Which.ALL)
