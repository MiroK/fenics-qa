from dolfin import *

# Let the tensor depend on domains
def tensor_conditional(predicate, tvalue, fvalue):
    '''
    tensor_conditional is a generalization of UFL.conditional. The latter allows
    the return values to be scalars only. In tensor_functional the true and 
    false value can have arbitrary (matching) shape.
    '''
    assert tvalue.ufl_shape == fvalue.ufl_shape

    rank = len(tvalue.ufl_shape)
    # Scalar case
    if rank == 0: return conditional(predicate, tvalue, fvalue)
    # Vector base case
    if rank == 1:
        dim = tvalue.ufl_shape[0]
        conds = [conditional(predicate, tvalue[i], fvalue[i]) for i in range(dim)]
        return as_vector(conds)
    # Higher order recursively
    else:
        dim = tvalue.ufl_shape[0]
        return as_vector([tensor_conditional(predicate, tvalue[i, Ellipsis], fvalue[i, Ellipsis])
                          for i in range(dim)])


i, j, k, l = indices(4)
delta = Identity(3)
C_1 = as_tensor(delta[i,j]*delta[k,l], (i, j, k, l))  # value for domain 1
C_2 = as_tensor(delta[i,k]*delta[j,l], (i, j, k, l))  # ... for domain 2
C_3 = as_tensor(delta[i,l]*delta[j,k], (i, j, k, l))  # ... for domain 2


mesh = UnitCubeMesh(4, 4, 6)
# The domains will be marked by DG0 function
S = FunctionSpace(mesh, 'DG', 0)
domains = interpolate(Expression('(x[2] > 0.5-DOLFIN_EPS) ? 1 : ((x[2] > 0.25-DOLFIN_EPS) ? 2 : 3)', degree=1), S)

C = tensor_conditional(lt(domains, 1.5),
                       C_1,
                       tensor_conditional(lt(domains, 2.5),
                                          C_2,
                                          C_3))

V = VectorFunctionSpace(mesh, 'P', 1)
du = TrialFunction(V)
del_u = TestFunction(V)
u = Function(V)

eps = sym(grad(u))
tau = as_tensor(C[j,i,k,l]*eps[k,l] , (j,i))

F_u = tau[j,i]*del_u[i].dx(j)*dx

Form = F_u 
Gain = derivative(Form, u, du)

A = assemble(Gain)
assert A.norm('linf') > 0
print A.norm('linf')
