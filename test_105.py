from dolfin import *
import numpy as np

mesh = UnitSquareMesh(10, 10)
gdim = mesh.geometry().dim()
V = VectorFunctionSpace(mesh, 'Quadrature', 2)
# Quad points
xq = V.dofmap().tabulate_all_coordinates(mesh).reshape((-1, gdim))
# You might want to remove the duplicates
xq0 = xq[V.sub(0).dofmap().dofs()]
