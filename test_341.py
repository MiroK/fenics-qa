from dolfin import *
from petsc4py import PETSc

mesh = UnitCubeMesh(3, 3, 3)

DG = FiniteElement('DG', mesh.ufl_cell(), 0)
Nedelec = FiniteElement('N1curl', mesh.ufl_cell(), 1)
element = MixedElement([DG, Nedelec])
V = FunctionSpace(mesh, element)

V_RT = FunctionSpace(mesh, 'RT', 1)

bc1 = DirichletBC(V_RT, Constant([0.0, 0.0, 0.0]), DomainBoundary())
bc2 = DirichletBC(V.sub(0), Constant(0.0), DomainBoundary())
bc3 = DirichletBC(V.sub(1), Constant([0.0, 0.0, 0.0]), DomainBoundary())

u = TrialFunction(V_RT)
q, eta = TestFunctions(V)

a = div(u)*q*dx + inner(u, curl(eta))*dx
L = q*dx
dummy = q*dx

A, _ = assemble_system(a, dummy, [bc1, bc2, bc3])
_, b = assemble_system(a, L, [bc1, bc2])

from block.iterative import CGN

AAinv = CGN(A, show=2)

x = AAinv*b
