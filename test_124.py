from dolfin import *
from numpy.linalg import eigvals

deg = 2

mesh = UnitCubeMesh(1, 1, 1)

# Where displacement, gradient and eigenvalues live
V = VectorFunctionSpace(mesh, 'CG', deg)
T = TensorFunctionSpace(mesh, 'DG', deg-1)
E = VectorFunctionSpace(mesh, 'DG', deg-1)

# Displacement field
v = interpolate(Expression(('x[0]', '2*x[1]', '3*x[2]'), degree=1), V)
# Gradient field
t = project(grad(v), T)
# Field of eigenvalues
e = Function(E)
e_values = e.vector().array()
# Buidl the eigenvalues from cell values
t_values = t.vector().array()
for cell in cells(mesh):
    cell_grads = t_values[T.dofmap().cell_dofs(cell.index())]
    cell_grads = cell_grads.reshape((9, T.element().space_dimension()/9)).T
    
    eigw_cell_dofs = E.dofmap().cell_dofs(cell.index())
    eigw_cell_dofs = eigw_cell_dofs.reshape((3, E.element().space_dimension()/3)).T

    for cell_grad, eigw_dofs in zip(cell_grads, eigw_cell_dofs):
        print eigvals(cell_grad.reshape((3, 3))).real
    print
# Finalize
e.vector()[:] = e_values

e0 = Constant((1, 2, 3))

e0 = interpolate(e0, E)
e0.vector().axpy(-1, e.vector())
print e0.vector().norm('linf')

print e0.vector().array()


# v = interpolate(Expression(('sin(pi*x[0])', 'sin(2*pi*x[1])', 'sin(3*pi*x[2])'), degree=1), V)
# e0 = Expression(('pi*cos(pi*x[0])', '2*pi*cos(2*pi*x[1])', '3*pi*cos(3*pi*x[2])'))
