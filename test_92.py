from dolfin import *

mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, 'CG', 1)
bc = DirichletBC(V, Constant(0), 'near(x[0], 0) || near(x[0], 1)')

gdim = mesh.geometry().dim()
d2v = dof_to_vertex_map(V)
for dof in bc.get_boundary_values():
    vertex = Vertex(mesh, d2v[dof])
    print [vertex.x(i) for i in range(gdim)]
