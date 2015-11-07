from dolfin import *

mesh = UnitSquareMesh(10, 10)
# Space for tensor components
S = FunctionSpace(mesh, 'DG', 0)
# Space for tensor
T = TensorFunctionSpace(mesh, 'DG', 0)
# Coordinates of dofs of S
dof_xs = S.dofmap().tabulate_all_coordinates(mesh).reshape((-1, 2))

c0, c1, c2, c3 = [Function(S)]*4
# Build the components
for i, dof_x in enumerate(dof_xs):
    x, y = dof_x
    c0.vector()[i] = x**2 + y**2 + 1
    c1.vector()[i] = x**2
    c2.vector()[i] = y**2
    c3.vector()[i] = 1

# Alternatively do the above by iterating over cells
dofmap = S.dofmap()
for cell in cells(mesh):
    dof = dofmap.cell_dofs(cell)[0]
    x, y = cell.midpoint().x(), cell.midpoint().y()
    c0.vector()[i] = x**2 + y**2 + 1
    c1.vector()[i] = x**2
    c2.vector()[i] = y**2
    c3.vector()[i] = 1

# Make the tensor
t = interpolate(Expression((('c0', 'c1'), ('c2', 'c3')),
                c0=c0, c1=c1, c2=c2, c3=c3),
                T)
