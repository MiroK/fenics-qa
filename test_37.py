from dolfin import *

mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)
M = PETScMatrix()
assemble(inner(u, v)*dx, M)

Mmat = M.mat()
print Mmat.getInfo()
print Mmat.sizes
print Mmat.size
print Mmat.local_size
print Mmat.block_size
print Mmat.block_sizes
print Mmat.owner_range
print Mmat.owner_ranges
print Mmat.assembled
print Mmat.symmetric
print Mmat.hermitian
print Mmat.structsym
