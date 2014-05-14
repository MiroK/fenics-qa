from dolfin import *

mesh = UnitSquareMesh(10, 10)
mesh_ = BoundaryMesh(mesh, 'exterior')

V = FunctionSpace(mesh, 'CG', 2)
V_ = FunctionSpace(mesh, 'CG', 2)

u_ = interpolate(Expression('sin(pi*x[0])*cos(pi*x[1])'), V_)

d = mesh.topology().dim()
bdr_to_full = mesh_.entity_map(d-1)

dofmap = V.dofmap()
u = Function(V)

for bdr_entity
