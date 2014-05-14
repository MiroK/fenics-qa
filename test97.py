from dolfin import *

mesh = UnitSquareMesh(2, 2)

V = VectorFunctionSpace(mesh, 'CG', 2)
W = FunctionSpace(mesh, 'CG', 1)
X = MixedFunctionSpace([V, W])

u = interpolate(Expression(('0', '0', '1')), X)

dofs_W = X.sub(1).dofmap().dofs()

print u.vector()[dofs_W].array()
