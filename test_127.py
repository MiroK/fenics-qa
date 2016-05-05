from dolfin import *

def tensor_conditional(predicate, tvalue, fvalue):
    assert tvalue.ufl_shape == fvalue.ufl_shape
    
    # Scalar
    rank = len(tvalue.ufl_shape)
    if rank == 0: 
        return conditional(predicate, tvalue, fvalue)
    # Vector
    elif rank == 1:
        dim = tvalue.ufl_shape[0]
        conds = [conditional(predicate, tvalue[i], fvalue[i]) for i in range(dim)]
    # Matrices
    elif rank == 2:
        nrows, ncols = tvalue.ufl_shape
        conds = [[conditional(predicate, tvalue[i, j], fvalue[i, j]) for j in range(ncols)]
                 for i in range(nrows)]
    # Generalize later
    else:
        raise ValueError
        
    return as_tensor(conds)

# ----------------------------------------------------------------------------

mesh = UnitSquareMesh(10, 10)

V = VectorFunctionSpace(mesh, 'CG', 1)
W = FunctionSpace(mesh, 'DG', 0)

w = interpolate(Expression('x[0]*x[0]+x[1]*x[1]'), W)

x, y = SpatialCoordinate(mesh)

try:
    b = conditional(gt(w, 1), as_vector((x, y)), as_vector((-x, -y)))
    f = project(b, V)
except Exception as e:
    print '\033[1;37;31m%s\033[0m' % e
    print 'Try with tensor_conditional'

b = tensor_conditional(gt(w, 1), as_vector((x, y)), as_vector((-x, -y)))
f = project(b, V)
plot(f)

V = FunctionSpace(mesh, 'CG', 1)
b = tensor_conditional(gt(w, 1),
                       as_tensor(((1, 0), (0, 2))),
                       Constant(1+1E-2)*as_tensor(((2, 0), (0, 1))))
f = project(det(b), V)
plot(f)

interactive()
