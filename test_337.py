from dolfin import *

mesh = UnitSquareMesh(20, 20)

V = VectorFunctionSpace(mesh, 'P', 2)
Q = FunctionSpace(mesh, 'P', 1)

PPoint  = 'near(x[0], 0.0) && near(x[1], 0.0)'
bcp = DirichletBC(Q, Constant(0.0), PPoint, 'pointwise')

px_ = Function(Q)
py_ = Function(Q)

p  = TrialFunction(Q)
v = TestFunction(V)
F5 = dot(as_vector([px_,py_]) - grad(p), v)*dx
a  = lhs(F5)
L  = rhs(F5)

A = assemble(a)
print A.size(0), A.size(1)
