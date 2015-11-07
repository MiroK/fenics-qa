from dolfin import *
import numpy as np

mesh = UnitSquareMesh(64, 64)

V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

f = Expression('(1+2*pi*pi)*cos(pi*x[0])*cos(pi*x[1])')
a = inner(grad(u), grad(v))*dx + inner(u, v)*dx
L = inner(f, v)*dx

uh = Function(V)
solve(a == L, uh)

# plot(uh, interactive=True)

dofmap = V.dofmap()
mesh.init(2, 0)
for cell in cells(mesh):
    cell_dofs = dofmap(cell.index())
