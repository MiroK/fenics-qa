from dolfin import *

mesh = UnitIntervalMesh(4)

V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

m = inner(u, v)*dx
m_1 = (1./CellSize(mesh))*inner(u, v)*dx(metadata={'quadrature_degree': 0})

M = assemble(m)
M_1 = assemble(m_1)

print M.array()
print M_1.array()

for cell in cells(mesh):
    print cell.volume()
