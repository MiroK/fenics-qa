from dolfin import *
parameters['reorder_dofs_serial'] = True

mesh = UnitIntervalMesh(4)

V1 = FunctionSpace(mesh, 'CG', 1)
V2 = FunctionSpace(mesh, 'DG', 0)
R = FunctionSpace(mesh, 'R', 0)
u = TrialFunction(V1)
X = MixedFunctionSpace([V2, R])
(v, c) = TestFunctions(X)

# this bilinear form does not involve the real test function c,
# so one row should be identically zero
a0 = u.dx(0) * v * dx
print assemble(a0).array()

# this bilinear form only involves the real test function c
# so the corresponding row should have a nonzero
a1 = u * c * Expression("1 - x[0]") * ds
print assemble(a1).array()

print 'V2'
print V2.dofmap().tabulate_all_coordinates(mesh)

print 'R'
print R.dofmap().tabulate_all_coordinates(mesh)

print 'X'
print X.dofmap().tabulate_all_coordinates(mesh)
print X.sub(0).dofmap().dofs()
print X.sub(1).dofmap().dofs()
