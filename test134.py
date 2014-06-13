from dolfin import *

mesh = UnitSquareMesh(3, 3)
# Mark boundaries of original mesh
bdry = FacetFunction('size_t', mesh, 0)
DomainBoundary().mark(bdry, 10)

# Function space for displacement u
V = VectorFunctionSpace(mesh, 'CG', 1)
# Shear displacement u
u = interpolate(Expression(('0.25*x[1]', '0')), V)

plot(bdry, title='Initial', interactive=True)

# Mote the mesh
mesh.move(u)
plot(bdry, title='Moved', interactive=True)

# Revert manually
x, y = mesh.coordinates()[:, 0], mesh.coordinates()[:, 1]
x -= 0.25*y
plot(bdry, title='And back', interactive=True)
