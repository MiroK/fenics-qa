from dolfin import *

dim = 2

mesh = UnitSquareMesh(5,5)

S = FunctionSpace(mesh, 'DG', 1)
V = VectorFunctionSpace(mesh, 'DG', 1)
W = TensorFunctionSpace(mesh, 'DG', 1, shape=(dim,dim))

U = interpolate(Expression(('x[0]', 'x[1]')), V)
M = Function(W)

for i in range(dim):
    for j in range(dim):
        print type(U[i].dx(j))
        ind = i*dim + j 
        assign(M.sub(ind), project(U[i].dx(j), S))

for i in range(dim):
    for j in range(dim):
        plot(M[i, j], title='[%d, %d]' % (i, j))
interactive()
