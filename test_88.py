from dolfin import *
import numpy as np

dt = 0.5
mesh = UnitSquareMesh(4, 4)
bmesh = RectangleMesh(Point(1, 1), Point(2, 2), 4, 4)
bmesh = BoundaryMesh(bmesh, 'exterior')
u = ALE.move(mesh, bmesh)

# Velocity as function
V = VectorFunctionSpace(mesh, 'CG', 1)
v = interpolate(u, V)
v.vector()[:] /= dt
bcf = DirichletBC(V, v, 'on_boundary')

# Velocity as expression
class Velocity(Expression):
    def __init__(self, u, dt):
        self.u = u
        self.dt = dt
    def eval_cell(self, value, x, cell):
        self.u.eval_cell(value, x, cell)
        value /= dt
    def value_shape(self):
        return (2, )

v = Velocity(u, dt)
bcg = DirichletBC(V, v, 'on_boundary')

# Check
f = bcf.get_boundary_values() 
g = bcg.get_boundary_values()
assert len(set(g.keys())-set(f.keys())) == 0
assert all(near(f[key], g[key]) for key in f)
