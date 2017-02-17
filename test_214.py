from matplotlib.tri import Triangulation, LinearTriInterpolator
from dolfin import *
import numpy as np


npoints = 20
# 'Measurement' points
x, y = np.tile(np.linspace(0, 5, 10), 2), np.repeat(np.linspace(0, 5, 10), 2)
# Delaunay triangulation for the points (considery meshpy, mshr?)
triang = Triangulation(x, y)
# Turn it into FEniCS mesh
mesh = Mesh()
mesh_vertices = np.c_[x, y]
mesh_cells = triang.triangles

editor = MeshEditor()
editor.open(mesh, 2, 2)
editor.init_vertices(len(mesh_vertices))
editor.init_cells(len(mesh_cells))

for vertex_index, v in enumerate(mesh_vertices): editor.add_vertex(vertex_index, v)

for cell_index, c in enumerate(mesh_cells): editor.add_cell(cell_index, *c)

editor.close()

# Now lets have data at the points
ux, uy = x**2+y**2, x*y

# First option to represent the data: Expression with callback to matplotlib
class TriInterp(Expression):
    def __init__(self, triang, ux, uy, **kwargs):
        self.Ix = LinearTriInterpolator(triang, ux)
        self.Iy = LinearTriInterpolator(triang, uy)

    def eval(self, value, x): 
        value[0] = self.Ix(*x)
        value[1] = self.Iy(*x)

    def value_shape(self): return (2, )

f = TriInterp(triang, ux, uy, degree=1, cell=mesh.ufl_cell())
L = inner(grad(f), grad(f))*dx(domain=mesh)
print assemble(L)

# For later use it will be faster to represent f as function
V = VectorFunctionSpace(mesh, 'CG', 1)
v = interpolate(f, V)

plot(v, interactive=True)
