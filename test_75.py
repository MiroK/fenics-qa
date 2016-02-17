from dolfin import *

mesh0 = UnitSquareMesh(1, 1)

mesh1 = RectangleMesh(Point(0, 0), Point(4, 4), 2, 2)

plot(mesh0)
plot(mesh1)
interactive()
