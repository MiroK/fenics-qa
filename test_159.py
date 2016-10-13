from dolfin import *
import ufl

mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

c0, c1 = Constant(1), Constant(2)

a = c0*inner(u, v)*dx + c0*c1*inner(grad(u), grad(v))*dx

for e in ufl.algorithms.iter_expressions(a):
    print e
    for t in ufl.algorithms.traverse_terminals(e):
        print '\t', t
