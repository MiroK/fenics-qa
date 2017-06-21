from dolfin import *
from mshr import *
import sys

domain = Sphere(Point(0, 0, 0), 2)
mesh = generate_mesh(domain, int(sys.argv[1])) 

V = VectorFunctionSpace(mesh, 'CG', 1)

# Projection of the normal vector on P2 space
u = TrialFunction(V)
v = TestFunction(V)
n = FacetNormal(mesh)

a = inner(u, v)*ds + Constant(0)*inner(u, v)*dx
L = inner(n, v)*ds

# Solve system
A = assemble(a)
b = assemble(L)
A.ident_zeros()

n = Function(V)
solve(A, n.vector(), b)

V = FunctionSpace(mesh, 'DG', 0)
u = TrialFunction(V)
v = TestFunction(V)

a = inner(u, v)*ds + Constant(0)*inner(u, v)*dx
L = inner(div(n), v)*ds

A = assemble(a)
b = assemble(L)
A.ident_zeros()

divn = Function(V)
solve(A, divn.vector(), b)

# File('x.pvd') << n
plot(divn)
interactive()

# n0 = Expression(('x[0]', 'x[1]', 'x[2]'), degree=1, cell=tetrahedron)
print assemble(divn*ds(domain=mesh))
