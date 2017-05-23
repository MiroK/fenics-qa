import sympy as sp
from dolfin import *

x, y, z = sp.symbols('x[0] x[1] x[2]')

u = x**2 + 2*y**2 + 3*z**2
ddu = sum(u.diff(xi, 2) for xi in (x, y, z))

u = Expression(sp.printing.ccode(u), degree=2, cell=tetrahedron)
ddu = Expression(sp.printing.ccode(ddu), degree=2, cell=tetrahedron)

DDU = div(grad(u))
e = DDU - ddu

mesh = UnitCubeMesh(4, 4, 4)
print assemble(inner(e, e)*dx(domain=mesh))



