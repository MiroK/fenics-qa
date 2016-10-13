from dolfin import *
from petsc4py import PETSc

mesh = UnitSquareMesh(100, 100)
V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

a = inner(grad(u), grad(v))*dx + inner(u, v)*dx
L = inner(Constant(1), v)*dx
A, b = assemble_system(a, L)

comm = mesh.mpi_comm().tompi4py()

Amat = as_backend_type(A).mat()
ksp = PETSc.KSP().create()
ksp.setOperators(Amat)
pc = ksp.getPC()
pc.setType(PETSc.PC.Type.HYPRE)

class MyPC(PETScUserPreconditioner):
    def __init__(self, pc):
        self.pc=pc
        PETScUserPreconditioner.__init__(self)

    def solve(self, x, b):
        print '.'
        self.pc.apply(as_backend_type(b).vec(),
                      as_backend_type(x).vec())

my_pc = MyPC(pc)
solver = PETScKrylovSolver('cg', my_pc)
solver.set_operator(A)
solver.parameters['monitor_convergence']=True
solver.parameters['relative_tolerance'] = 1E-13
solver.parameters['absolute_tolerance'] = 1E-13

uh = Function(V)
x = uh.vector()
niters = solver.solve(x, b)
