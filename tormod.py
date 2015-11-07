from dolfin import *

mesh = UnitTriangleMesh()
V = FunctionSpace(mesh, 'CG', 1)
Q = FunctionSpace(mesh, 'CG', 3)

u = TrialFunction(V)
v = TestFunction(Q)
a = inner(u, v)*dx

A = assemble(a)
n_rows, n_cols = A.size(0), A.size(1)
assert n_cols == V.dim() and n_rows == Q.dim()

# Element matrix
# Test
n0 = a.arguments()[0].function_space().dolfin_element().space_dimension()
# Trial
n1 = a.arguments()[1].function_space().dolfin_element().space_dimension()

assert n_rows == n0 and n_cols == n1

mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, 'CG', 1)
Q = FunctionSpace(mesh, 'CG', 3)

u = TrialFunction(V)
v = TestFunction(Q)
a = inner(u, v)*dx

assert n0 == a.arguments()[0].function_space().dolfin_element().space_dimension()
assert n1 == a.arguments()[1].function_space().dolfin_element().space_dimension()
