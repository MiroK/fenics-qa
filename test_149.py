from dolfin import *
from itertools import product

mesh = BoxMesh(Point(-1, -1, -1), Point(1, 1, 1), 2, 2, 3)

domains = CellFunction('size_t', mesh, 0)
qs = [CompiledSubDomain('x[0] %s= 0 && x[1] %s= 0' % (s, t))
      for s, t in product(('>', '<'), ('>', '<'))]

for i, q in enumerate(qs, 1): q.mark(domains, i)

submeshes = [SubMesh(mesh, domains, i) for i in range(1, 5)]

[submesh.move(Expression(('a', 'b', '0'), a=a, b=b))
for (submesh, (a, b)) in zip(submeshes, product((0.2, -0.2), (0.2, -0.2)))]

for i, mesh in enumerate(submeshes, 1):
    f = XDMFFile(mpi_comm_world(), "q%d_mesh.xdmf" % i)
    f << mesh
