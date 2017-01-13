from dolfin import *

mesh = UnitSquareMesh(10, 10)
V = VectorElement("Lagrange", mesh.ufl_cell(), 1, dim=2)

ME = FunctionSpace(mesh, V*V*V)
for i in range(ME.num_sub_spaces()): print ME.sub(i).num_sub_spaces()
print

ME = FunctionSpace(mesh, MixedElement([V, V, V]))
for i in range(ME.num_sub_spaces()): print ME.sub(i).num_sub_spaces()
