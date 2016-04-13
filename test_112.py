from dolfin import *


mesh = UnitSquareMesh(10 ,10)
V = FunctionSpace(mesh, 'CG', 1 )
v = TestFunction(V)
L = Constant(0.0)*v*dx
b = assemble(L)

list_of_points = Point(0.25, 0.25), Point(0.5, 0.5), Point(0.75, 0.75)
amp = 2.

for pt in list_of_points:
    print 'b before', b.norm('linf')
    source_add = PointSource(V, Point(pt), amp)
    source_add.apply(b)
    print 'b with add', b.norm('linf')
    # Use this for some calculation.
    # remove the PointSource so I don't need to reassemble b.
    source_remove = PointSource(V, Point(pt), -amp)
    source_remove.apply(b)
    print 'b after', b.norm('linf')
