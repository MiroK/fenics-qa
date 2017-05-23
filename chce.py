# coding: utf-8
from dolfin import *
mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, 'CG', 1)
dofmap = V.dofmap()
print dofmap.dofs()
print dofmap.tabulate_local_to_global_dofs()
