from dolfin import *

mesh = UnitSquareMesh(2, 2)

f = Expression('sin(pi*x[0])*x[1]')
V = FunctionSpace(mesh, 'DG', 2)

f_proj = project(f, V)

F = f_proj.vector().array()
X = V.dofmap().tabulate_all_coordinates(mesh)
X.resize((V.dim(), 2))

print 'dof index | dof coordinate |  dof value'
for i, (x, v) in enumerate(zip(X, F)):
    print i, x, v
