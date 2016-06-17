from dolfin import *

mesh = RectangleMesh(Point(0,0), Point(1,2), 1,2)

sd = CompiledSubDomain("x[1] <= 1.0")
mf = MeshFunction("size_t", mesh, 2)
mf.set_all(0)
sd.mark(mf, 1)

top = SubMesh(mesh, mf, 0)
bottom = SubMesh(mesh, mf, 1)
down = Expression(["0.0", "-1.0"])

bottom.move(down)
bottomfile = XDMFFile(mpi_comm_world(), "bottom.xdmf")
topfile = XDMFFile(mpi_comm_world(), "top.xdmf")

bottomfile << bottom
topfile << top
