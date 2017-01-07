from dolfin import *

mesh = UnitSquareMesh(10,10)
P1 = FiniteElement('P', triangle, 1)
P = P1*P1
V = FunctionSpace(mesh, P)

# view
try:
    view = V.sub(0)
    Function(view)
except RuntimeError:
    print 'Using proper FunctionSpace'
    V1 = V.sub(0).collapse()
    Function(V1)

dim = V1.dim()
assert view.dim() == dim

# Illustrate that V.sub(0) is related to mixed function space by showing that
# some dofs (local) exceed dim
assert all(V1.dofmap().dofs() < dim)

print any(view.dofmap().dofs() > dim)
