from dolfin import *

mesh = UnitSquareMesh(2, 2)
V = VectorFunctionSpace(mesh, 'CG', 1)
v = TestFunction(V)
b = assemble(inner(Constant((0, 0)), v)*dx)

x = Point(0.25, 0.25)
p0 = PointSource(V.sub(0), x, 1)
p1 = PointSource(V.sub(1), x, 2)

p0.apply(b)
print b.array()
p1.apply(b)
print b.array()
