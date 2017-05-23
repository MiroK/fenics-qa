from dolfin import *
from block import block_mat, block_vec
import numpy as np

mesh = UnitSquareMesh(10, 10)
V = VectorFunctionSpace(mesh, 'CG', 1)

u, v = TrialFunction(V), TestFunction(V)
dofs_x = V.sub(0).dofmap().dofs()
dofs_y = V.sub(1).dofmap().dofs()

M = assemble(inner(div(u), div(v))*dx)

# With cbc.block
Q = FunctionSpace(mesh, 'CG', 1)
p, q = TrialFunction(Q), TestFunction(Q)
block00 = assemble(inner(p.dx(0), q.dx(0))*dx)
block01 = assemble(inner(p.dx(1), q.dx(0))*dx)
block10 = assemble(inner(p.dx(0), q.dx(1))*dx)
block11 = assemble(inner(p.dx(1), q.dx(1))*dx)
Mb = block_mat([[block00, block01],
                [block10, block11]])

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
