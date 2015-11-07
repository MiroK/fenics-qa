from dolfin import *

p = Point(1, 2, 3)
x, y, z = p[:3]

print x, y, z
