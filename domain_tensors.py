from dolfin import *
import ufl


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


def nested_if(predicates, statements):
    '''
    if predicates[0]:
        statements[0]
    else if predicates[1]:
        statements[1]
    else:
        statements[-1]
    '''
    assert (len(predicates)+1) == len(statements)
    assert len(set([s.ufl_shape for s in statements])) == 1
    if len(predicates) == 1:
        return conditional(predicates[0], statements[0], statements[1])
    else:
        return conditional(predicates[0],
                           statements[0],
                           nested_if(predicates[1:], statements[1:]))


def switch(expr, cases, statements):
    '''
    This is a swith statement. It reads as follows

    swich(expr):
    case: cases[0]   # expr == cases[0]
        statements[0]
    case: cases[1]   # expr == cases[1]
        statements[1]
    default:
        statements[-1]
    '''
    cases = [eq(expr, c) for c in cases]
    return nested_if(cases, statements)

# -----------------------------------------------------------------------------

if __name__ == '__main__':
    from itertools import product

    mesh = RectangleMesh(Point(-1, -1), Point(1, 1), 10, 10)
    V = FunctionSpace(mesh, 'DG', 0)

    f = interpolate(Expression('x[0]'), V)
    g = project(nested_if([lt(f, 0.1), lt(f, 0.8)],
                          [Constant(1), Constant(2), Constant(3)]), V)
    h = project(switch(g**2, (1, 4), map(Constant, (10, 20, 0))), V)

    V = FunctionSpace(mesh, 'CG', 1)
    u, v = TrialFunction(V), TestFunction(V)
    C = tensor_conditional(h > 2, Constant((1, 2)), Constant((1, -3)))

    a = tr(outer(C, C))*inner(u, v)*dx
    A = assemble(a)
