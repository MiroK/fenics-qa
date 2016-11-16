from dolfin import *

mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)
a = inner(grad(u), grad(v))*dx
A = PETScMatrix(); assemble(a, A)
A_h = PETScMatrix(); assemble(a, A_h)

bc = DirichletBC(V, Constant(0), 'on_boundary')
bc_h = DirichletBC(V, Constant(1), 'on_boundary')
foo = bc_h.homogenize()      # Now bc_h acts like bc
print type(foo)              # Not DirichletBC but none!

bc.apply(A)
print '|A|', A.norm('linf')
bc_h.apply(A_h)
print '|A_h|', A_h.norm('linf')
# Show they are same
A -= A_h
print '|A-A_h|', A.norm('linf')
