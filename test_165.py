from mpi4py import MPI as pyMPI
import numpy as np
from dolfin import BoxMesh, Point

angle = 46.
L = 4.6

mesh = BoxMesh(Point(-L, -L, -L), Point(L, L, L), 10, 10, 10)
mesh.rotate(angle, 2, Point(0., 0., 0.))

x = mesh.coordinates()

comm = mesh.mpi_comm().tompi4py()
# Local
local_min = np.min(x, axis=0)
local_max = np.max(x, axis=0)
# Alloc global
global_min = np.tile(local_min, (comm.size, 1))
global_max = np.tile(local_max, (comm.size, 1))
# Communicate
comm.Allgather([local_min, pyMPI.DOUBLE], [global_min, pyMPI.DOUBLE])
comm.Allgather([local_max, pyMPI.DOUBLE], [global_max, pyMPI.DOUBLE])
# 'Reduce'
bounding_min = np.min(global_min, axis=0)
bounding_max = np.max(global_max, axis=0)

# plot(mesh, interactive=True)

if comm.rank == 0:
    print bounding_min, bounding_max
    
    from dolfin import between

    if between(angle, (0, 90)):
        from dolfin import cos, pi, sqrt

        angle = np.deg2rad(angle)
        A = sqrt(2)*L*cos(pi/4-angle)

        min0 = np.array([-A, -A, -L])
        max0 = np.array([A, A, L])

        assert np.linalg.norm(min0-bounding_min, np.inf) < 1E-13
        assert np.linalg.norm(max0-bounding_max, np.inf) < 1E-13
