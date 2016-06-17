from dolfin import *

mesh = UnitSquareMesh(10, 10)
V = VectorFunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

e = Constant((1, 2))
A = assemble(inner(dot(grad(u), e), v)*dx)
B = A.copy()
A_mat, B_mat = as_backend_type(A).mat(), as_backend_type(B).mat()

x = interpolate(Expression(('sin(x[0])', 'cos(x[1])'), degree=2), V).vector()

#--- Transpose
# By action
y1 = x.copy()
A.transpmult(x, y1)

# With explicit matrix B = A'
A_mat.transpose(B_mat)

# We have non-sym matrix so sanity check
A -= B
e = A.norm('linf')
info('|A-B|=%g' % e)

# Make sure they are the same
y2 = x.copy()
B.mult(x, y2)

y2.axpy(-1, y1)
e = y2.norm('linf')
info('|y1-y2|=%g' % e)

#--- Power
# By action
for i in range(2): A.mult(x, y1)

# With explicit matrix
C = A.copy(); C_mat = as_backend_type(C).mat()
C_mat.matMult(C_mat)  # C = C**2
C.mult(x, y2)

# Make sure they are the same
y2.axpy(-1, y1)
e = y2.norm('linf')
info('|y1-y2|=%g' % e)
