from dolfin import *
from mshr import *
import numpy as np

domain = Rectangle(Point(-1, -1), Point(1, 1))
mesh = generate_mesh(domain, 20)

class Metric(Expression):
    def __init__(self, mesh, **kwargs):
        self.mesh = mesh

    def eval_cell(self, values, x, cell):
        x = Cell(self.mesh, cell.index).get_vertex_coordinates()
        x0, x1, x2 = x.reshape(3, 2)
        # [x, = x0*(1-s-t) + x1*s + x2*t = x0 + [x1-x0 x2-x0][s,
        #  y]                                                 t]
        F = np.c_[x1-x0, x2-x0]   # d(x, y)/d(s, t)
        Finv = np.linalg.inv(F)
        values[:] = (Finv.T.dot(Finv)).flatten()

    def value_shape(self): return (2, 2)

f = Metric(mesh, degree=0)
M = TensorFunctionSpace(mesh, 'DG', 0)
fh = interpolate(f, M)

for cell in cells(mesh):
    v0 = cell.volume()
    v = 0.5/sqrt(np.linalg.det(fh(cell.midpoint()).reshape(2, 2)))
    assert near(v, v0)
