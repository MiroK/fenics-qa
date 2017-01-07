from dolfin import *
from petsc4py import PETSc
from mpi4py import MPI as pyMPI
import numpy as np

mesh = RectangleMesh(Point(1, 1), Point(2, 2), 40, 40)
nsub = 4

element = MixedElement([FiniteElement('Lagrange', mesh.ufl_cell(), i)
                        for i in range(1, 1+nsub)])

W = FunctionSpace(mesh, element)

f = tuple(['x[0]*x[0]+x[1]*x[1]+%d' % i for i in range(1, 1+nsub)])
f = Expression(f, degree=2)
f = interpolate(f, W)
vec = f.vector()

comm = mesh.mpi_comm().tompi4py()
# Find min by PETSc
if parameters['linear_algebra_backend'] == 'PETSc':
    petsc_vec = as_backend_type(vec).vec()

    for i in range(nsub):
        Wi_dofs = W.sub(i).dofmap().dofs()
        Wi_is = PETSc.IS().createGeneral(Wi_dofs, comm)

        sub_vec = petsc_vec.getSubVector(Wi_is)
        minimum = sub_vec.min()[-1]
        info('PETSc, sub %d, Min is %g' % (i, minimum))

# Alternatively w/out PETSc
first, last = W.dofmap().ownership_range()
for i in range(nsub):
    Wi_dofs = W.sub(i).dofmap().dofs()
    Wi_dofs = Wi_dofs[np.logical_and(Wi_dofs > first, Wi_dofs <= last)] - first

    values = vec.get_local()
    # Owned & of subspace i
    values = values[Wi_dofs]

    local_min = np.min(values)
    # Reduce over processes
    minimum = comm.allreduce(local_min, op=pyMPI.MIN)
    info('noPETSc, sub %d, Min is %g' % (i, minimum))
