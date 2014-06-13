from dolfin import *
import numpy as np

if has_linear_algebra_backend('uBLAS'):
    parameters['linear_algebra_backend'] = 'uBLAS'

    A = uBLASDenseMatrix(2, 2)
    cols0 = np.array([0, 1], dtype='uintp')
    vals0 = np.array([3., 2.])
    cols1 = np.array([0, 1], dtype='uintp')
    vals1 = np.array([-2., -1.])
    A.setrow(0, cols0, vals0)
    A.setrow(1, cols1, vals1)
    A.apply('insert')
    print 'A =', A.array()

    b = uBLASVector(2)
    b_ = np.array([3., -1.])
    b.set_local(b_)
    print 'b =', b.array()

    x = uBLASVector(2)
    solve(A, x, b)
    print 'x =', x.array()
