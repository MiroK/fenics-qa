from dolfin import *
from mshr import *

domain = Circle(Point(0, 0), 1)
mesh = generate_mesh(domain, 80)

mesh = BoundaryMesh(mesh, 'exterior')
print 'PI is', sum(cell.volume() for cell in cells(mesh))/2
