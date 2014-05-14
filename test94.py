from dolfin import *

mesh = UnitSquareMesh(10, 10)

V = VectorFunctionSpace(mesh, 'CG', 1)

u = interpolate(Expression(('1 + x[0]', '2 + x[1]')), V)
u_v = u.vector()

dof_coordinates = V.dofmap().tabulate_all_coordinates(mesh)
n = V.dim()
d = mesh.geometry().dim()
dof_coordinates.resize((n, d))

x_dofs = V.sub(0).dofmap().dofs()
y_dofs = V.sub(1).dofmap().dofs()

for x_dof, y_dof in zip(x_dofs, y_dofs):
  print dof_coordinates[x_dof], dof_coordinates[y_dof], u_v[x_dof], u_v[y_dof]
