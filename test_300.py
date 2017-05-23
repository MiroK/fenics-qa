from dolfin import *
from block import block_mat, block_vec
import numpy as np

mesh = UnitSquareMesh(10, 10)
V = VectorFunctionSpace(mesh, 'CG', 1)

v = TestFunction(V)
vx, vy = v[0], v[1]  # In the assembly these act like (vx, 0), (0, vy)
dofs_x = V.sub(0).dofmap().dofs()
dofs_y = V.sub(1).dofmap().dofs()

f = Constant(1)

# (vx, 0) - so second comp is 0
b = assemble(inner(vx, f)*dx)
b_array = b.array()
print np.linalg.norm(b_array[dofs_y])
# (0, vy) - so first comp is 0
b = assemble(inner(vy, f)*dx)
b_array = b.array()
print np.linalg.norm(b_array[dofs_x])

# Consider the mass matrix
u = TrialFunction(V)
M = assemble(inner(u, v)*dx)

# With cbc.block
Q = FunctionSpace(mesh, 'CG', 1)
p, q = TrialFunction(Q), TestFunction(Q)
block = assemble(inner(p, q)*dx)
Mb = block_mat([[block, 0],
                [0, block]])

# The equivalence: let's see how the two operators act on representation of the same
# function. For cbc block the dofs of subspaces are serialized.
x_array = np.random.rand(V.dim())
x = Function(V).vector(); x.set_local(x_array); x.apply('insert')
# Reorder for cbc.block
xx_array = x_array[dofs_x]
xy_array = x_array[dofs_y]
xx = Function(Q).vector(); xx.set_local(xx_array); xx.apply('insert')
xy = Function(Q).vector(); xy.set_local(xy_array); xy.apply('insert')
xb = block_vec([xx, xy])

# Action
y = x.copy()
M.mult(x, y)
yx, yy = y.array()[dofs_x], y.array()[dofs_y]

yb = Mb*xb
ybx, yby = yb[0].array(), yb[1].array()

print np.linalg.norm(yx-ybx), np.linalg.norm(yy-yby)
