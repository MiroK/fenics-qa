from dolfin import *

rect0 = Rectangle(-1, -1, 0, 1)
rect1 = Rectangle(0, -1, 1, 1)
domain = rect0 + rect1
mesh = Mesh(domain, 20)

plot(mesh)
interactive()
