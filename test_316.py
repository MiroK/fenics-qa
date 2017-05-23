from dolfin import *
import sympy as sp

# Compute the symbolic expression for eigenvalues by sympy
T = sp.Matrix(2, 2, lambda i, j: sp.Symbol('T[%d, %d]' % (i, j), real=True))
eig_expr = T.eigenvects()   # ((v, multiplicity, [w])...)

eigv = [e[0] for e in eig_expr]
eigw = [e[-1][0] for e in eig_expr]

eigv_expr = map(str, eigv)
eigw_expr = [[str(e[0]), str(e[1])] for e in eigw]

# UFL operator for eigenvalues of 2x2 matrix, a pair of scalars
def eigv(T): return map(eval, eigv_expr)

# UFL operator for eigenvectors of 2x2 matrix, a pair of vectors
def eigw(T): return [as_vector(map(eval, vec)) for vec in eigw_expr]

n = 50
mesh = UnitSquareMesh(n, n)

V = TensorFunctionSpace(mesh, 'CG', 1)
x = SpatialCoordinate(mesh)

t = as_matrix(((x[0]+1, x[0]+x[1]), (x[0]+x[1], x[1]+1)))
t = project(t, V)

v0, v1 = eigv(t)
w0, w1 = eigw(t)

# Eigenvalues okay
print sqrt(assemble(inner(det(t)-v0*v1, det(t)-v0*v1)*dx))
print sqrt(assemble(inner(tr(t)-(v0+v1), tr(t)-(v0+v1))*dx))

# Orthogonality of eigenvectors
print sqrt(assemble(inner(Constant(0) - dot(w0, w1), Constant(0) - dot(w0, w1))*dx))

# Decomposition
e = dot(t, w0)
e0 = v0*w0
print sqrt(assemble(inner(e-e0, e-e0)*dx))

e = dot(t, w1)
e0 = v1*w1
print sqrt(assemble(inner(e-e0, e-e0)*dx))
