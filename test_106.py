from dolfin import *

n = 10

mesh = UnitSquareMesh(n, n)
V = FunctionSpace(mesh, 'CG', 1)
psi = Function(V)

def u0_boundary(x,on_boundary):
    return on_boundary
bc = DirichletBC(V,Constant(0.0),u0_boundary)

u = TrialFunction(V)
v = TestFunction(V)

a = inner(Constant(300.)*grad(u), grad(v))*dx + inner(u, v)*dx
m = inner(u, v)*dx
L = Constant(0.)*v*dx

A, M, _ = PETScMatrix(), PETScMatrix(), PETScVector()

assemble_system(a,L,bc,A_tensor=A,b_tensor=_)
assemble_system(m,L,A_tensor=M,b_tensor=_)

print A.sparray()

eigensolver = SLEPcEigenSolver(A,M)
eigensolver.parameters['spectrum'] = 'target real'
eigensolver.parameters['tolerance'] = 1.e-14
eigensolver.parameters['problem_type'] = 'gen_hermitian'
eigensolver.parameters['spectral_transform'] = 'shift-and-invert'
eigensolver.parameters['spectral_shift'] = 5924.

eigensolver.solve(5)

assert eigensolver.get_number_converged() > 0

for i in range(eigensolver.get_number_converged()):
    r, c, rx, cx = eigensolver.get_eigenpair(i)
    print "%d:  %g" % (0, r)
    plot(Function(V, rx), title='%d' % i)
interactive()
