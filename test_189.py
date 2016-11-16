from dolfin import *

mesh = UnitSquareMesh(10, 10)
V = FiniteElement('Lagrange', mesh.ufl_cell(), 1)
V = MixedElement([V, V])
V = FunctionSpace(mesh, V)

u = Function(V.sub(0).collapse())
u.assign(Constant(4))
as_backend_type(u.vector()).vec().
