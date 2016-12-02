from dolfin import *

mesh = UnitSquareMesh(10, 10)
V = VectorFunctionSpace(mesh, 'CG', 2)
u = interpolate(Expression(('x[0]', 'x[0]*x[1]'), degree=2), V)

F = Identity(2) + grad(u)
C = F*F.T

# Eigenvalues are roots of characteristic polynomial
# e**2 - tr(A)*e + det(A)
def eig_plus(A): return (tr(A) + sqrt(tr(A)**2-4*det(A)))/2

def eig_minus(A): return (tr(A) - sqrt(tr(A)**2-4*det(A)))/2

# Check
S = FunctionSpace(mesh, 'CG', 2)

f0 = project(tr(C), S)
f = project(eig_plus(C)+eig_minus(C), S)
print (f0.vector()-f.vector()).norm('l2')

f0 = project(det(C), S)
f = project(eig_plus(C)*eig_minus(C), S)
print (f0.vector()-f.vector()).norm('l2')
