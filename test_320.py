from dolfin import *

mesh = UnitSquareMesh(10, 10)

M = MixedElement([
    MixedElement([FiniteElement('BDM', mesh.ufl_cell(), 1),
                  FiniteElement('BDM', mesh.ufl_cell(), 1)]),
                  VectorElement('DG', mesh.ufl_cell(), 0),
                  FiniteElement('DG', mesh.ufl_cell(), 0) ])
W = FunctionSpace(mesh, M)
test = TestFunctions(W)

#tau = as_vector((test[0], test[1]))
#v = test[-2]
#eta = test[-1]

tau, v, eta = TrialFunctions(W)
#tau = as_vector((tau[0], tau[1]))

w_n = Function(W)
w_n.vector()[:] = 0.

sigma_n, u_n, gamma_n = w_n.split()
#sigma_n = as_vector((sigma_n[0], sigma_n[1]))
#sigma_n = as_vector((sigma0, sigma1))
#u_n = Function(VectorFunctionSpace(mesh, 'DG', 0))
#sigma_n = Function(VectorFunctionSpace(mesh, 'BDM', 1))

Meta = as_matrix([[0., eta],[-eta, 0.]])

b_below = dot(v, div(sigma_n))*dx #+ 
b_below += inner(Meta, sigma_n)*dx

assemble(b_below)
