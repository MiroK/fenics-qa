from dolfin import *
from sympy import Matrix
mesh = UnitSquareMesh(10, 10)
V = VectorFunctionSpace(mesh, 'CG', 2)
u = interpolate(Expression(('x[0]', 'x[0]*x[1]'), degree=2), V)

F = Identity(2) + grad(u)
C = F*F.T

# Eigenvalues are roots of characteristic polynomial
# e**2 - tr(A)*e + det(A) = 0
def eig_plus(A): 
    return (tr(A) + sqrt(tr(A)**2-4*det(A)))/2
def eig_minus(A): 
    return (tr(A) - sqrt(tr(A)**2-4*det(A)))/2
def eig_vecmat(A):
    lambdap = eig_plus(A)
    lambdan = eig_minus(A)
    a = A[0,0]
    b = A[0,1]
    c = A[1,0]
    d = A[1,1]

    return conditional(ne(c, 0),
                       as_matrix(((lambdap-d,lambdan-d), (c,c))),
                       conditional(ne(b, 0),
                                   as_matrix(((b,b), (lambdap-a,lambdan-a))),
                                   Identity(2)))
# Check
S = FunctionSpace(mesh, 'CG', 2)
T = TensorFunctionSpace(mesh, 'CG', 2)
f0 = project(eig_vecmat(C),T) 
plot(f0[0,0],interactive=True)
