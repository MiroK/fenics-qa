from dolfin import *
import sys

n = int(sys.argv[1])
mesh = UnitCubeMesh(n, n, n)
cell = mesh.ufl_cell()

V = FiniteElement('RT', cell, 2)
Q = VectorElement('Lagrange', cell, 1)
W = MixedElement([V, Q])

W = FunctionSpace(mesh, W)

u, p = TrialFunctions(W)
v, q = TestFunctions(W)

x, y, z = SpatialCoordinate(mesh)
f = x*sin(pi*y)+z
a = inner(div(u), div(v))*dx + inner(v, curl(p))*dx + inner(u, curl(q))*dx
L = inner(f, div(v))*dx

bc = DirichletBC(W.sub(0), Constant((0, 0, 0)), 'on_boundary')

A, b = assemble_system(a, L, bc)

wh = Function(W)
solve(A, wh.vector(), b)

uh, ph = wh.split()

print sqrt(assemble(inner(uh, uh)*dx))
print sqrt(assemble(inner(ph, ph)*dx))

plot(uh)
interactive()
