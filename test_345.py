from fenics import *
import numpy as np
from numpy import linalg

mesh = UnitCubeMesh(10, 10, 10)
RT = FunctionSpace(mesh, 'RT', 1)
DG = FunctionSpace(mesh, 'DG', 0)
Nedelec = FunctionSpace(mesh, 'N1curl', 1)

u = TrialFunction(RT)
q = TestFunction(DG)
eta = TestFunction(Nedelec)

bc_RT = DirichletBC(RT, Constant([0.0, 0.0, 0.0]), DomainBoundary())
phi = Expression('x[0]+x[1]', degree=1)

a = div(u)*q*dx
L = phi*q*dx
a_ = inner(u, curl(eta))*dx
dummy = inner(Constant([0,0,0]), eta)*dx

A, b = assemble_system(a, L)
A_, _ = assemble_system(a_, dummy)

AA = np.vstack((A.array(), A_.array())) # merge two blocks of matrix(stiffness matrix)
zero = np.zeros(A_.array().shape[0])
bb = np.hstack((b.array(), zero)) # merge two blocks of vector(source term)

uu = linalg.lstsq(AA, bb)[0] # use least square algorithm in numpy

u = Function(RT)
u.vector().set_local(uu)

plot(u, title='solution u')

