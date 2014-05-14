from dolfin import *
import sys

#define mesh and function space
mesh = CircleMesh(Point(0,0),5,0.7)
V = FunctionSpace(mesh, 'Lagrange', 1)

#build essential boundary conditions
def u0_boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V,Constant(0.0) , u0_boundary)

#define functions
u = TrialFunction(V)
v = TestFunction(V)

Pot = Expression('0.0')

#define problem
a = (inner(grad(u), grad(v)) \
     + Pot*u*v)*dx
m = u*v*dx

A = PETScMatrix()
M = PETScMatrix()
_ = PETScVector()
L = Constant(0.)*v*dx

assemble_system(a, L, bc, A_tensor=A, b_tensor=_)

if sys.argv[1] == 'use':
    assemble_system(m, L, bc, A_tensor=M, b_tensor=_)
else:
    assemble_system(m, L, A_tensor=M, b_tensor=_)

#create eigensolver
eigensolver = SLEPcEigenSolver(A,M)
eigensolver.parameters['spectrum'] = 'smallest magnitude'
eigensolver.parameters['spectral_shift'] = 1.0
eigensolver.parameters['tolerance'] = 1.e-15

#solve for eigenvalues
eigensolver.solve()

u = Function(V)
for i in range(eigensolver.get_number_converged()):
    #extract next eigenpair
    r, c, rx, cx = eigensolver.get_eigenpair(i)
    if not near(r, 1, 1E-10):
        print 'eigenvalue:', r, 'eigenfunction max', rx.max()

        #assign eigenvector to function
        u.vector()[:] = rx
        plot(u,interactive=True, title='real %d' % i)
