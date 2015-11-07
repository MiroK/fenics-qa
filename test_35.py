from dolfin import *
import numpy as np
np.set_printoptions(precision=2)

mesh = UnitSquareMesh(2, 2)
V = FunctionSpace(mesh, 'CG', 1)
W = MixedFunctionSpace([V, V])

u0, u1 = TrialFunctions(W)
v0, v1 = TestFunction(W)

M = assemble(inner(u0, v1)*dx)
print 'M'
print M.array()

c = interpolate(Constant((1, 1)), W)
m_lumped = assemble(action(inner(u0, v1)*dx, c))
M_lumped = assemble(Constant(0)*inner(u0, v1)*dx)
M_lumped.zero()
M_lumped.set_diagonal(m_lumped)
print 'M lumped'
print M_lumped.array()

# W_dofmap = W.dofmap()
# dof_x = W.dofmap().tabulate_all_coordinates(mesh).reshape((-1, 2))
# 
# print 'W[0] dofs'
# for dof in W.sub(0).dofmap().dofs():
#     print dof, dof_x[dof]
# 
# print 'W[1] dofs'
# for dof in W.sub(1).dofmap().dofs():
#     print dof, dof_x[dof]
# 
